#!/usr/bin/env python3
"""葉緑体ゲノム gb → ビューワー用の注釈データ (data/genomes.js) を生成する。

各ゲノムについて:
  - 自己 blastn で IR(逆位反復)を検出し、LSC が先頭(0)に来るよう座標を回転
  - CDS/tRNA/rRNA を遺伝子ブロックとして抽出（相同遺伝子をそろえる正規化キー付き）
  - 機能カテゴリ分類、トランススプライシング遺伝子(psaA 等)は最大エキソンで代表

出力 data/genomes.js は `window.GENOME_INDEX` と `window.GENOME_DATA` を定義する
（file:// で直接開けるよう JSON ではなく JS で埋め込む）。

依存: Biopython, NCBI BLAST+ (makeblastdb, blastn)
使い方: python3 generate_data.py
"""
import json
import re
import subprocess
import tempfile
import time
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
SRC = HERE / "gb"                  # gb はこのフォルダで一括管理（無ければアクセッションから自動取得）
SRC.mkdir(exist_ok=True)
OUT = HERE / "data"
OUT.mkdir(exist_ok=True)

# シンテニー（相同ヒット）として保存・表示する最小値。ここがサイトのスライダー下限にもなる。
BLAST_MIN_PID = 30
BLAST_MIN_COV = 25

# (label, 表示名, gbファイル名, アクセッション)
GENOMES = [
    ("Nicotiana_benthamiana",     "Nicotiana benthamiana",     "nicotiana_benthamiana.gb",     "cultivar LAB"),
    ("Chlamydomonas_reinhardtii", "Chlamydomonas reinhardtii", "chlamydomonas_reinhardtii.gb", "NC_005353.1"),
    ("Manihot_esculenta",         "Manihot esculenta (cassava)", "manihot_esculenta.gb",       "NC_010433.1"),
    ("Sorghum_bicolor",           "Sorghum bicolor (sorghum)",   "sorghum_bicolor.gb",         "NC_008602.1"),
    ("Cryptomeria_japonica",      "Cryptomeria japonica (sugi)", "cryptomeria_japonica.gb",    "NC_010548.1"),
    ("Marchantia_polymorpha",     "Marchantia polymorpha (liverwort)", "marchantia_polymorpha.gb", "NC_037507.1"),
    ("Arabidopsis_thaliana",      "Arabidopsis thaliana",            "arabidopsis_thaliana.gb",      "NC_000932.1"),
    ("Ceratopteris_thalictroides", "Ceratopteris thalictroides (water fern)", "ceratopteris_thalictroides.gb", "NC_062137.1"),
    ("Oryza_sativa",              "Oryza sativa (rice)",             "oryza_sativa.gb",              "NC_001320.1"),
    ("Cyanidioschyzon_merolae",   "Cyanidioschyzon merolae (red alga)", "cyanidioschyzon_merolae.gb", "NC_004799.1"),
    ("Zea_mays",                  "Zea mays (maize)",                "zea_mays.gb",                  "NC_001666.2"),
    ("Phaeodactylum_tricornutum", "Phaeodactylum tricornutum (diatom)", "phaeodactylum_tricornutum.gb", "NC_008588.1"),
]


def ensure_gb(path, acc):
    """gb が無ければ NCBI からアクセッションで取得（accession が NCBI 形式のときのみ）。"""
    if path.exists() and path.stat().st_size > 0:
        return
    if not re.match(r"^[A-Z]{1,2}_?\d", acc):      # 'cultivar LAB' のような非アクセッションは取得不可
        raise SystemExit(f"missing gb and not a fetchable accession: {path.name} ({acc})")
    print(f"  fetching {acc} -> gb/{path.name}")
    for _ in range(5):
        with open(path, "w") as fh:
            subprocess.run(["efetch", "-db", "nuccore", "-id", acc, "-format", "gb"], stdout=fh)
        if path.stat().st_size > 0 and path.read_text(errors="ignore").startswith("LOCUS"):
            return
        time.sleep(3)
    raise SystemExit(f"failed to fetch {acc}")

# ---- 遺伝子名の正規化（相同遺伝子をそろえる）----
RRNA_KEY = {
    "16s": "rrn_16S", "rrn16": "rrn_16S", "rrns": "rrn_16S",
    "23s": "rrn_23S", "rrn23": "rrn_23S", "rrnl": "rrn_23S",
    "5s": "rrn_5S", "rrn5": "rrn_5S",
    "4.5s": "rrn_4.5S", "rrn4.5": "rrn_4.5S",
    "7s": "rrn_7S", "rrn7": "rrn_7S",
    "3s": "rrn_3S", "rrn3": "rrn_3S",
}
SPLIT_GENE = {"rpob1": "rpob", "rpob2": "rpob", "rpoc1a": "rpoc1", "rpoc1b": "rpoc1"}


def normalize(gene, product, ftype):
    g = (gene or "").strip().lower()
    p = (product or "").strip().lower()
    if ftype == "tRNA":
        m = re.match(r"^trn([a-z]+)", g)
        if m:
            return "trn" + m.group(1)
        m = re.search(r"trna-([a-z]{3})", p)
        return "trn_" + m.group(1) if m else (g or p)
    if ftype == "rRNA":
        for token in (g, re.sub(r"\s*ribosomal.*", "", p)):
            t = token.replace(" ", "")
            if t in RRNA_KEY:
                return RRNA_KEY[t]
        return g or p
    key = g.replace("_", "").replace("-", "").replace(" ", "") or p
    return SPLIT_GENE.get(key, key)


def func_cat(key):
    if key.startswith("trn"): return "tRNA"
    if key.startswith("rrn"): return "rRNA"
    if key.startswith("ndh"): return "ndh"
    if key.startswith("psa"): return "PSI"            # 光化学系I
    if key.startswith("psb"): return "PSII"           # 光化学系II
    if key.startswith("pet"): return "cytb6f"         # シトクロム b6f
    if key.startswith("atp") or key == "rbcl": return "photosynthesis"   # ATP合成酵素 / RuBisCO
    if key.startswith(("rps", "rpl")): return "ribosomal"
    if key.startswith("rpo"): return "RNApol"
    if key.startswith("chl"): return "chl"
    return "other"


def revcomp(s):
    return s.translate(str.maketrans('ACGTUNacgtun', 'TGCAANtgcaan'))[::-1]


def ir_duplicate(genes, regions):
    """IR は逆位反復で2コピーある。注釈は片方(通常IRa)のみなので、もう片方へ
    遺伝子をミラー（位置反転＋鎖反転）して複製する。key に __ir2 を付けて区別。"""
    irs = {r["label"]: r for r in regions if r["label"] in ("IRa", "IRb")}
    if "IRa" not in irs or "IRb" not in irs:
        return []
    ra, rb = irs["IRa"], irs["IRb"]
    K = ra["start"] + rb["end"]                    # mirror(p) = K - p  (IRa[a0+i] <-> IRb[b1-i])

    def in_region(g, r):
        m = (g["start"] + g["end"]) / 2
        return r["start"] <= m < r["end"]

    cnt_a = sum(in_region(g, ra) for g in genes)
    cnt_b = sum(in_region(g, rb) for g in genes)
    src = ra if cnt_a >= cnt_b else rb
    if max(cnt_a, cnt_b) == 0:
        return []
    out = []
    for g in genes:
        if g["cat"] == "intergenic" or not in_region(g, src):
            continue
        out.append({**g, "key": g["key"] + "__ir2", "start": K - g["end"],
                    "end": K - g["start"], "strand": -g["strand"]})
    return out


def compute_intergenic(genes, length, label, min_len=100):
    """遺伝子ブロックの隙間（遺伝子間領域 IGR）を遺伝子と同じ形式で返す。座標は回転後。"""
    prot = [g for g in genes]
    merged = []
    for s, e in sorted((g["start"], g["end"]) for g in prot):
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    gaps, cur = [], 0
    for s, e in merged:
        if s - cur >= min_len:
            gaps.append((cur, s))
        cur = max(cur, e)
    if length - cur >= min_len:
        gaps.append((cur, length))

    def left_name(a):
        c = [g for g in prot if g["end"] <= a]
        return max(c, key=lambda g: g["end"])["name"] if c else "start"

    def right_name(b):
        c = [g for g in prot if g["start"] >= b]
        return min(c, key=lambda g: g["start"])["name"] if c else "end"

    out = []
    for i, (s, e) in enumerate(gaps):
        out.append({"key": f"{label}__ig{i}", "name": f"IGR {left_name(s)}–{right_name(e)}",
                    "start": s, "end": e, "strand": 1, "cat": "intergenic"})
    return out


def gene_nt(seq_rot, gene):
    s = seq_rot[gene['start']:gene['end']]
    return revcomp(s) if gene['strand'] < 0 else s


def feature_extent(f, max_span=10000):
    s, e = int(f.location.start), int(f.location.end)
    if e - s > max_span and len(f.location.parts) > 1:
        p = max(f.location.parts, key=lambda q: int(q.end) - int(q.start))
        return int(p.start), int(p.end)
    return s, e


# ---- IR 検出と四分割構造 ----
def detect_regions(seq, length, tmp, ir_min_bp=5000):
    fa = tmp / "g.fasta"
    with open(fa, "w") as fh:
        fh.write(">g\n")
        for i in range(0, length, 70):
            fh.write(seq[i:i + 70] + "\n")
    subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl", "-out", str(tmp / "g")],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = subprocess.run(
        ["blastn", "-task", "blastn", "-query", str(fa), "-db", str(tmp / "g"),
         "-evalue", "1e-10", "-outfmt", "6 qstart qend sstart send length"],
        capture_output=True, text=True, check=True).stdout
    irs = []
    for line in out.strip().splitlines():
        qs, qe, ss, se, ln = map(int, line.split("\t"))
        if ln >= ir_min_bp and ss > se:
            irs.append((min(qs, qe), max(qs, qe), min(ss, se), max(ss, se)))
    if not irs:
        return [(0, length, "LSC")]
    irs.sort(key=lambda t: t[1] - t[0], reverse=True)
    a0, a1, b0, b1 = irs[0]
    irx, iry = sorted([(a0, a1), (b0, b1)])
    gap1, gap2 = iry[0] - irx[1], (length - iry[1]) + irx[0]
    g1, g2 = ("LSC", "SSC") if gap1 >= gap2 else ("SSC", "LSC")
    segs = []
    if irx[0] > 0:
        segs.append((0, irx[0], g2))
    segs += [(irx[0], irx[1], "IRa"), (irx[1], iry[0], g1), (iry[0], iry[1], "IRb")]
    if iry[1] < length:
        segs.append((iry[1], length, g2))
    return segs


def lsc_start(segs, length):
    """LSC が先頭・連続になる回転原点を返す。LSC が原点をまたぐ場合はその周回開始位置。"""
    lsc = [(s, e) for s, e, lab in segs if lab == "LSC"]
    if not lsc:
        return 0
    wraps = any(e == length for s, e in lsc) and any(s == 0 for s, e in lsc)
    if wraps:                                   # 原点をまたぐ→末尾側ピースの開始を原点に
        return next(s for s, e in lsc if e == length)
    return lsc[0][0]


def rot_interval(a, b, off, L):
    if off == 0:
        return (a, b)
    if a < off < b:
        return None
    ra, rb = (a - off) % L, (b - off) % L
    return (ra, rb if rb != 0 else L)


def rotate_segments(segs, off, L):
    out = sorted((r[0], r[1], lab) for s, e, lab in segs if (r := rot_interval(s, e, off, L)))
    merged = []
    for s, e, lab in out:
        if merged and merged[-1][2] == lab and abs(merged[-1][1] - s) < 1:
            merged[-1] = (merged[-1][0], e, lab)
        else:
            merged.append((s, e, lab))
    return [{"start": int(s), "end": int(e), "label": lab} for s, e, lab in merged]


def build_genome(label, gb_path, tmp):
    rec = SeqIO.read(gb_path, "genbank")
    seq = str(rec.seq)
    L = len(seq)
    regions_raw = detect_regions(seq, L, tmp)
    off = lsc_start(regions_raw, L)
    regions = rotate_segments(regions_raw, off, L)

    genes, seen = [], set()
    for f in rec.features:
        if f.type not in ("CDS", "tRNA", "rRNA"):
            continue
        gene = f.qualifiers.get("gene", [None])[0]
        product = f.qualifiers.get("product", [None])[0]
        key = normalize(gene, product, f.type)
        if key in seen:
            continue
        seen.add(key)
        s, e = feature_extent(f)
        rs, re_ = (s - off) % L, (e - off) % L
        if re_ < rs:                                   # 回転の切れ目をまたぐ→中点で近似
            mid = ((s + e) // 2 - off) % L; w = e - s
            rs, re_ = mid - w / 2, mid + w / 2
        genes.append({
            "key": key,
            "name": gene or product or key,
            "start": int(min(rs, re_)),
            "end": int(max(rs, re_)),
            "strand": int(f.location.strand or 1),
            "cat": func_cat(key),
        })
    genes.sort(key=lambda g: g["start"])
    genes += ir_duplicate(genes, regions)          # IRb にも遺伝子をミラー複製
    genes.sort(key=lambda g: g["start"])
    genes += compute_intergenic(genes, L, label)   # 遺伝子間領域(IGR)も追加（cat=intergenic）
    genes.sort(key=lambda g: g["start"])
    seq_rot = (seq[off:] + seq[:off]).upper()      # 配列も同じ回転を適用（pos 0 = LSC開始）
    return {"length": L, "regions": regions, "genes": genes}, seq_rot


def main():
    index, data, seqs = [], {}, {}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for label, disp, fname, acc in GENOMES:
            gb = SRC / fname
            ensure_gb(gb, acc)                      # 無ければアクセッションから自動取得
            g, seq_rot = build_genome(label, gb, tmp)
            g["label"] = label
            g["display"] = disp
            g["accession"] = acc
            data[label] = g
            seqs[label] = seq_rot
            index.append({"label": label, "display": disp, "length": g["length"],
                          "accession": acc, "nGenes": len(g["genes"])})
            print(f"{disp:34} {g['length']:>7,} bp  genes={len(g['genes']):>3}  "
                  f"regions={'/'.join(r['label'] for r in g['regions'])}")

        # --- per-gene BLAST (all ordered pairs) for the strict synteny mode ---
        import itertools
        gene_fa = {}
        for label in data:
            fa = tmp / f"{label}.genes.fasta"
            with open(fa, "w") as fh:
                for gene in data[label]["genes"]:
                    nt = gene_nt(seqs[label], gene)
                    if len(nt) >= 20:
                        fh.write(f">{gene['key']}\n{nt}\n")
            gene_fa[label] = fa
        from collections import defaultdict
        cen = {lab: {g["key"]: (g["start"] + g["end"]) / 2 / data[lab]["length"] for g in data[lab]["genes"]} for lab in data}
        blast = {}
        print("\nper-gene blastn (strict synteny):")
        for A, B in itertools.permutations(list(data), 2):
            db = tmp / f"{B}.genedb"
            subprocess.run(["makeblastdb", "-in", str(gene_fa[B]), "-dbtype", "nucl", "-out", str(db)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out_b = subprocess.run(
                ["blastn", "-task", "blastn", "-query", str(gene_fa[A]), "-db", str(db),
                 "-evalue", "1e-3", "-max_target_seqs", "10",
                 "-outfmt", "6 qseqid sseqid pident length qlen bitscore"],
                capture_output=True, text=True).stdout
            hits = defaultdict(list)
            for line in out_b.strip().splitlines():
                q, s, pid, length, qlen, bits = line.split("\t")
                pid, length, qlen, bits = float(pid), int(length), int(qlen), float(bits)
                cov = round(length / qlen * 100, 1) if qlen else 0
                if pid >= BLAST_MIN_PID and cov >= BLAST_MIN_COV:
                    hits[q].append((bits, pid, cov, s))
            out_list = []
            for q, lst in hits.items():
                topb = max(h[0] for h in lst)
                cands = [h for h in lst if h[0] >= topb * 0.99]   # bitscore がほぼ同率の候補
                fq = cen[A].get(q, 0.0)
                # 同率なら「位置が近い＝同じ側のIRコピー」を選ぶ（交差リボンを防ぐ）
                bits, pid, cov, s = min(cands, key=lambda h: (abs(cen[B].get(h[3], 0.0) - fq), -h[0]))
                out_list.append({"q": q, "s": s, "pid": round(pid, 1), "cov": cov})
            blast[f"{A}::{B}"] = out_list
            print(f"  {A:26}-> {B:26} {len(out_list):>3} gene hits")

    out = OUT / "genomes.js"
    with open(out, "w") as fh:
        fh.write("// Auto-generated by generate_data.py — do not edit by hand.\n")
        fh.write("window.GENOME_INDEX = " + json.dumps(index, ensure_ascii=False) + ";\n")
        fh.write("window.GENOME_DATA = " + json.dumps(data, ensure_ascii=False) + ";\n")

    seq_out = OUT / "sequences.js"               # 配列は別ファイル（やや大きいため）
    with open(seq_out, "w") as fh:
        fh.write("// Auto-generated by generate_data.py — rotated genome sequences.\n")
        fh.write("window.GENOME_SEQ = " + json.dumps(seqs, ensure_ascii=False) + ";\n")

    blast_out = OUT / "blast.js"                  # 厳密版シンテニー用の遺伝子ペアBLAST
    with open(blast_out, "w") as fh:
        fh.write("// Auto-generated by generate_data.py — per-gene blastn best hits per ordered pair.\n")
        fh.write(f"window.GENOME_BLAST_META = {{\"minPid\": {BLAST_MIN_PID}, \"minCov\": {BLAST_MIN_COV}}};\n")
        fh.write("window.GENOME_BLAST = " + json.dumps(blast, ensure_ascii=False) + ";\n")
    print(f"\nwrote {out} ({out.stat().st_size/1024:.0f} KB), "
          f"{seq_out} ({seq_out.stat().st_size/1024:.0f} KB), "
          f"{blast_out} ({blast_out.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
