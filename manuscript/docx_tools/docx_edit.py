"""Tracked, purple-coloured edits with explanatory comments for an unpacked .docx (after merge_runs.py).

Doc(unpacked_dir)
  .paras                      list of <w:p> elements in document order
  .text(i)                    visible text of paragraph i (ignores deleted text)
  .find(substr)               indices of paragraphs containing substr
  .replace(i, old, new, why)  tracked deletion of `old` (must be a contiguous substring of the paragraph's
                              visible text) + tracked purple insertion of `new`; comment `why` anchored to it
  .insert_after(i, text, why, style=None)  new tracked paragraph after paragraph i (purple), comment `why`
  .delete_para(i, why)        tracked deletion of the whole paragraph
  .comment(i, why, old=None)  comment only (anchored to `old` within paragraph i, or to the whole paragraph)
  .save()
"""
import copy
import sys
from datetime import datetime, timezone
from pathlib import Path
from lxml import etree

SKILL = Path("/Users/kizie/Library/Application Support/Claude/local-agent-mode-sessions/skills-plugin/8879dd4d-1049-447a-89bb-74aba2ceb3dd/36e6cf99-bd73-4f25-b6e1-acaca24f5810/skills/docx/scripts")
sys.path.insert(0, str(SKILL))
from comment import add_comment  # noqa: E402

WNS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
XMLNS = "http://www.w3.org/XML/1998/namespace"
W = "{%s}" % WNS
AUTHOR = "Claude (AI assistant)"
INITIALS = "AI"
PURPLE = "7030A0"
DATE = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _el(tag, **attrs):
    e = etree.Element(W + tag)
    for k, v in attrs.items():
        e.set(W + k, str(v))
    return e


class Doc:
    def __init__(self, unpacked):
        self.dir = Path(unpacked)
        self.path = self.dir / "word" / "document.xml"
        self.tree = etree.parse(str(self.path))
        self.body = self.tree.getroot().find(W + "body")
        self.paras = list(self.body.iter(W + "p"))
        self._rev = 9000
        self.log = []
        self.orig = [self.text(i) for i in range(len(self.paras))]   # text before any edit; anchors are matched here
        self.dry = False          # True: only check that every anchor / old text exists, collect failures in self.failures
        self.failures = []

    # ---------- reading ----------
    def _runs(self, p):
        """text runs of a paragraph that are not inside a deletion"""
        out = []
        for r in p.iter(W + "r"):
            if r.getparent().tag == W + "del" or r.find(W + "t") is None:
                continue
            out.append(r)
        return out

    def text(self, i):
        return "".join(t.text or "" for r in self._runs(self.paras[i]) for t in r.findall(W + "t"))

    @staticmethod
    def _norm(z):
        return (z.replace("\xa0", " ").replace("\u2009", " ").replace("\u202f", " ").replace("\u2019", "'").replace("\u2018", "'")
                 .replace("\u201c", '"').replace("\u201d", '"'))

    def find(self, substr):
        substr = self._norm(substr)
        return [i for i in range(len(self.paras)) if substr in self._norm(self.orig[i])]

    def one(self, substr):
        hits = self.find(substr)
        if len(hits) != 1:
            if self.dry and not hits:
                self.failures.append(f"anchor not found: {substr[:70]!r}")
                return 0
            raise ValueError(f"{len(hits)} paragraphs contain {substr[:60]!r}")
        return hits[0]

    # ---------- primitives ----------
    def _nid(self):
        self._rev += 1
        return self._rev

    def _split_run(self, r, k):
        """split run r at character offset k of its text (text = all <w:t> of the run, in order; tabs, breaks and
        page-break markers have length 0 and go to the right part); returns (left, right)"""
        kids = [e for e in r if e.tag != W + "rPr"]
        total = sum(len(e.text or "") for e in kids if e.tag == W + "t")
        if k <= 0:
            return None, r
        if k >= total:
            return r, None
        r2 = copy.deepcopy(r)
        for e in [e for e in r2 if e.tag != W + "rPr"]:
            r2.remove(e)
        pos, moved = 0, False
        for e in kids:
            if moved:
                r2.append(e)                          # moves the element
                continue
            L = len(e.text or "") if e.tag == W + "t" else 0
            if e.tag == W + "t" and pos < k < pos + L:
                t2 = copy.deepcopy(e)
                e.text = (e.text or "")[: k - pos]
                t2.text = (t2.text or "")[k - pos:]
                for q in (e, t2):
                    q.set("{%s}space" % XMLNS, "preserve")
                r2.append(t2)
                moved = True
            elif pos + L == k and L > 0:
                moved = True                          # split exactly after this text node
            pos += L
        r.addnext(r2)
        return r, r2

    def _locate(self, p, old):
        runs = self._runs(p)
        pos, s = [], ""
        for r in runs:
            txt = "".join(t.text or "" for t in r.findall(W + "t"))
            pos.append((len(s), len(s) + len(txt)))
            s += txt
        norm = lambda z: z.replace("\xa0", " ").replace("\u2009", " ").replace("\u202f", " ").replace("\u2019", "'").replace("\u2018", "'").replace("\u201c", '"').replace("\u201d", '"')
        s = norm(s); old = norm(old)
        k = s.find(old)
        if k < 0:
            import difflib
            sm = difflib.SequenceMatcher(None, s, old); m = sm.find_longest_match(0, len(s), 0, len(old))
            raise ValueError(f"text not found in paragraph: {old[:70]!r}\n   longest common part ends at: ...{old[max(0, m.b + m.size - 40):m.b + m.size]!r} | paragraph continues: {s[m.a + m.size:m.a + m.size + 60]!r} | old continues: {old[m.b + m.size:m.b + m.size + 60]!r}")
        if s.find(old, k + 1) >= 0:
            raise ValueError(f"text occurs twice in paragraph: {old[:70]!r}")
        a, b = k, k + len(old)
        # split boundary runs so the match covers whole runs
        sel = []
        for r, (ra, rb) in zip(runs, pos):
            if rb <= a or ra >= b:
                continue
            if ra < a:
                _, r = self._split_run(r, a - ra)
                ra = a
            if rb > b:
                r, _ = self._split_run(r, b - ra)
            sel.append(r)
        return sel

    def _purple_run(self, like, text):
        r = _el("r")
        rpr = copy.deepcopy(like.find(W + "rPr")) if like is not None and like.find(W + "rPr") is not None else _el("rPr")
        for tag in ("color", "highlight", "strike", "dstrike", "u"):
            for e in rpr.findall(W + tag):
                rpr.remove(e)
        c = _el("color", val=PURPLE)
        # schema order: rStyle, rFonts, b, bCs, i, iCs, caps, smallCaps, strike, dstrike, outline, shadow, emboss,
        # imprint, noProof, snapToGrid, vanish, webHidden, color, ... -> insert colour before spacing/w/kern/position/sz
        after = {W + t for t in ("spacing", "w", "kern", "position", "sz", "szCs", "highlight", "u", "effect", "bdr", "shd",
                                 "fitText", "vertAlign", "rtl", "cs", "em", "lang", "eastAsianLayout", "specVanish", "oMath")}
        idx = len(rpr)
        for j, e in enumerate(rpr):
            if e.tag in after:
                idx = j
                break
        rpr.insert(idx, c)
        r.append(rpr)
        t = _el("t")
        t.text = text
        t.set("{%s}space" % XMLNS, "preserve")
        r.append(t)
        return r

    def _comment_markers(self, cid):
        cs = _el("commentRangeStart", id=cid)
        ce = _el("commentRangeEnd", id=cid)
        rr = _el("r")
        rpr = _el("rPr")
        rpr.append(_el("rStyle", val="CommentReference"))
        rr.append(rpr)
        rr.append(_el("commentReference", id=cid))
        return cs, ce, rr

    def _add_comment(self, why):
        cid, _, _ = add_comment(self.dir, why, author=AUTHOR, initials=INITIALS)
        return cid

    # ---------- public edits ----------
    def replace(self, i, old, new, why=None):
        p = self.paras[i]
        if self.dry:
            try:
                self._locate(p, old)
            except ValueError as e:
                self.failures.append(f"para {i}: {e}")
            return
        sel = self._locate(p, old)
        d = _el("del", id=self._nid(), author=AUTHOR, date=DATE)
        sel[0].addprevious(d)
        for r in sel:
            d.append(r)                       # moves the run into the deletion
            for t in r.findall(W + "t"):
                t.tag = W + "delText"
                t.set("{%s}space" % XMLNS, "preserve")
        ins = None
        if new:
            ins = _el("ins", id=self._nid(), author=AUTHOR, date=DATE)
            ins.append(self._purple_run(sel[0], new))
            d.addnext(ins)
        if why:
            cid = self._add_comment(why)
            cs, ce, rr = self._comment_markers(cid)
            d.addprevious(cs)
            last = ins if ins is not None else d
            last.addnext(ce)
            ce.addnext(rr)
        self.log.append(("replace", i, old[:50], new[:50]))

    def comment(self, i, why, old=None):
        p = self.paras[i]
        if self.dry:
            if old:
                try:
                    self._locate(p, old)
                except ValueError as e:
                    self.failures.append(f"para {i}: {e}")
            return
        cid = self._add_comment(why)
        cs, ce, rr = self._comment_markers(cid)
        if old:
            sel = self._locate(p, old)
            sel[0].addprevious(cs)
            sel[-1].addnext(ce)
        else:
            runs = self._runs(p)
            if not runs:
                p.append(cs)
                p.append(ce)
            else:
                runs[0].addprevious(cs)
                runs[-1].addnext(ce)
        ce.addnext(rr)
        self.log.append(("comment", i, (old or "")[:50], why[:40]))

    def insert_after(self, i, text, why=None, like_run=None, keep_style=True):
        ref = self.paras[i]
        p = _el("p")
        ppr = copy.deepcopy(ref.find(W + "pPr")) if (keep_style and ref.find(W + "pPr") is not None) else _el("pPr")
        for e in ppr.findall(W + "rPr"):
            ppr.remove(e)
        rpr = _el("rPr")
        rpr.append(_el("ins", id=self._nid(), author=AUTHOR, date=DATE))
        ppr.append(rpr)                       # tracked paragraph mark
        p.append(ppr)
        ins = _el("ins", id=self._nid(), author=AUTHOR, date=DATE)
        like = like_run if like_run is not None else (self._runs(ref)[0] if self._runs(ref) else None)
        ins.append(self._purple_run(like, text))
        p.append(ins)
        if why:
            cid = self._add_comment(why)
            cs, ce, rr = self._comment_markers(cid)
            ins.addprevious(cs)
            ins.addnext(ce)
            ce.addnext(rr)
        ref.addnext(p)
        self.paras.insert(i + 1, p)
        self.orig.insert(i + 1, "")
        self.log.append(("insert_after", i, text[:50], ""))
        return i + 1

    def delete_para(self, i, why=None):
        p = self.paras[i]
        if self.dry:
            return
        runs = self._runs(p)
        if runs:
            d = _el("del", id=self._nid(), author=AUTHOR, date=DATE)
            runs[0].addprevious(d)
            for r in runs:
                d.append(r)
                for t in r.findall(W + "t"):
                    t.tag = W + "delText"
                    t.set("{%s}space" % XMLNS, "preserve")
        ppr = p.find(W + "pPr")
        if ppr is None:
            ppr = _el("pPr")
            p.insert(0, ppr)
        rpr = ppr.find(W + "rPr")
        if rpr is None:
            rpr = _el("rPr")
            ppr.append(rpr)
        rpr.insert(0, _el("del", id=self._nid(), author=AUTHOR, date=DATE))
        if why:
            cid = self._add_comment(why)
            cs, ce, rr = self._comment_markers(cid)
            first = p.find(W + "del")
            if first is not None:
                first.addprevious(cs)
                first.addnext(ce)
                ce.addnext(rr)
        self.log.append(("delete_para", i, self.text(i)[:50], ""))

    def save(self):
        self.tree.write(str(self.path), xml_declaration=True, encoding="UTF-8", standalone=True)
