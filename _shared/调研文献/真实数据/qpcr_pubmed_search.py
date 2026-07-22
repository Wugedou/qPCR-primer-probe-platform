#!/usr/bin/env python3
"""
qPCR 8病原体 × 3轮 PubMed 文献检索脚本
使用 Bio.Entrez (esearch + efetch abstract)
输出结构化 JSON: qpcr_search_raw.json
"""

import json
import time
import os
import re
from datetime import datetime

from Bio import Entrez

# NCBI 要求填写邮箱
Entrez.email = "zhangxiaoguang@example.com"
Entrez.tool = "qpcr-lit-search"
Entrez.do_dtd = False  # 禁用 DTD 验证，避免网络错误

# ========== 8 病原体配置 ==========
PATHOGENS = [
    {
        "name": "SARS-CoV-2",
        "genes": "N/RdRp",
        "queries": {
            1: '(SARS-CoV-2 OR "severe acute respiratory syndrome coronavirus 2") AND qPCR AND ("standard curve" OR LOD OR "limit of detection") AND copies',
            2: '(SARS-CoV-2 OR "severe acute respiratory syndrome coronavirus 2") AND qPCR AND (primer AND probe AND sequence) AND (N OR RdRp)',
            3: '(SARS-CoV-2 OR "severe acute respiratory syndrome coronavirus 2") AND qPCR AND (specificity OR "cross-reactivity" OR matrix) AND (N OR RdRp)',
        }
    },
    {
        "name": "Influenza A",
        "genes": "M/HA",
        "queries": {
            1: '("Influenza A" OR "Influenza A virus") AND qPCR AND ("standard curve" OR LOD OR "limit of detection") AND copies',
            2: '("Influenza A" OR "Influenza A virus") AND qPCR AND (primer AND probe AND sequence) AND (M OR HA OR hemagglutinin)',
            3: '("Influenza A" OR "Influenza A virus") AND qPCR AND (specificity OR "cross-reactivity" OR matrix) AND (M OR HA)',
        }
    },
    {
        "name": "Dengue virus",
        "genes": "NS1/3'UTR",
        "queries": {
            1: '("Dengue virus" OR DENV) AND qPCR AND ("standard curve" OR LOD OR "limit of detection") AND copies',
            2: '("Dengue virus" OR DENV) AND qPCR AND (primer AND probe AND sequence) AND (NS1 OR 3UTR OR 3\'UTR)',
            3: '("Dengue virus" OR DENV) AND qPCR AND (specificity OR "cross-reactivity" OR matrix) AND (NS1 OR 3UTR)',
        }
    },
    {
        "name": "HBV",
        "genes": "X/S",
        "queries": {
            1: '(HBV OR "hepatitis B virus") AND qPCR AND ("standard curve" OR LOD OR "limit of detection") AND copies',
            2: '(HBV OR "hepatitis B virus") AND qPCR AND (primer AND probe AND sequence) AND (X OR S OR surface OR polymerase)',
            3: '(HBV OR "hepatitis B virus") AND qPCR AND (specificity OR "cross-reactivity" OR matrix)',
        }
    },
    {
        "name": "M. tuberculosis",
        "genes": "IS6110",
        "queries": {
            1: '("Mycobacterium tuberculosis" OR MTB) AND (IS6110 OR IS6110) AND qPCR AND ("standard curve" OR LOD OR "limit of detection") AND copies',
            2: '("Mycobacterium tuberculosis" OR MTB) AND (IS6110) AND qPCR AND (primer AND probe AND sequence)',
            3: '("Mycobacterium tuberculosis" OR MTB) AND (IS6110) AND qPCR AND (specificity OR "cross-reactivity" OR matrix)',
        }
    },
    {
        "name": "HIV-1",
        "genes": "gag/LTR",
        "queries": {
            1: '(HIV-1 OR "human immunodeficiency virus type 1") AND qPCR AND ("standard curve" OR LOD OR "limit of detection") AND copies',
            2: '(HIV-1 OR "human immunodeficiency virus type 1") AND qPCR AND (primer AND probe AND sequence) AND (gag OR LTR OR "long terminal repeat")',
            3: '(HIV-1 OR "human immunodeficiency virus type 1") AND qPCR AND (specificity OR "cross-reactivity" OR matrix) AND (gag OR LTR)',
        }
    },
    {
        "name": "ZIKV",
        "genes": "E/NS5",
        "queries": {
            1: '(ZIKV OR "Zika virus") AND qPCR AND ("standard curve" OR LOD OR "limit of detection") AND copies',
            2: '(ZIKV OR "Zika virus") AND qPCR AND (primer AND probe AND sequence) AND (E OR NS5 OR envelope)',
            3: '(ZIKV OR "Zika virus") AND qPCR AND (specificity OR "cross-reactivity" OR matrix)',
        }
    },
    {
        "name": "E. coli O157",
        "genes": "rfbE/stx2",
        "queries": {
            1: '("Escherichia coli O157" OR "E. coli O157" OR EHEC O157) AND qPCR AND ("standard curve" OR LOD OR "limit of detection") AND copies',
            2: '("Escherichia coli O157" OR "E. coli O157" OR EHEC O157) AND qPCR AND (primer AND probe AND sequence) AND (rfbE OR stx2 OR "Shiga toxin")',
            3: '("Escherichia coli O157" OR "E. coli O157" OR EHEC O157) AND qPCR AND (specificity OR "cross-reactivity" OR matrix)',
        }
    },
]

OUTPUT_FILE = os.path.join(
    os.path.dirname(__file__),
    "qpcr_search_raw.json"
)

MAX_RESULTS_PER_SEARCH = 30  # 每轮最多返回 30 篇
SLEEP_BETWEEN_REQUESTS = 0.5  # NCBI 速率限制（无 API key 约 3/sec）


def safe_esearch(query, retmax=MAX_RESULTS_PER_SEARCH):
    """带重试的 esearch 调用"""
    for attempt in range(3):
        try:
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=retmax,
                sort="relevance",
                usehistory="n",
            )
            results = Entrez.read(handle)
            handle.close()
            return results
        except Exception as e:
            print(f"  ⚠ esearch 失败 (attempt {attempt+1}/3): {e}")
            time.sleep(2 ** attempt)
    return None


def safe_efetch(pmid_list):
    """批量 efetch abstract，返回 dict {pmid: {title, abstract, journal, year}}"""
    if not pmid_list:
        return {}
    result = {}
    records = []
    for attempt in range(3):
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid_list,
                rettype="xml",
                retmode="xml",
            )
            records = Entrez.read(handle)
            handle.close()
            break
        except Exception as e:
            print(f"  ⚠ efetch 失败 (attempt {attempt+1}/3): {e}")
            time.sleep(2 ** attempt)
            records = []
    
    if not records:
        return {}
    
    # records 是 DictionaryElement，含 'PubmedArticle' 和 'PubmedBookArticle'
    articles = list(records.get('PubmedArticle', []))
    articles += list(records.get('PubmedBookArticle', []))
    
    for article in articles:
        try:
            medline = article.get("MedlineCitation", article)
            article_data = medline.get("Article", {})
            pmid = str(medline.get("PMID", ""))
            
            title = str(article_data.get("ArticleTitle", ""))
            
            abstract_obj = article_data.get("Abstract", {})
            abstract_texts = abstract_obj.get("AbstractText", [])
            if isinstance(abstract_texts, list):
                abstract = " ".join([str(a) for a in abstract_texts])
            else:
                abstract = str(abstract_texts)
            
            journal_obj = article_data.get("Journal", {})
            journal_name = str(journal_obj.get("Title", ""))
            
            # 提取年份
            year = None
            pub_date = journal_obj.get("JournalIssue", {}).get("PubDate", {})
            year_str = pub_date.get("Year", "")
            if year_str:
                try:
                    year = int(year_str)
                except:
                    pass
            if not year:
                # 尝试从 MedlineDate 提取
                medline_date = pub_date.get("MedlineDate", "")
                match = re.search(r"(\d{4})", medline_date)
                if match:
                    year = int(match.group(1))
            
            result[pmid] = {
                "title": title,
                "abstract": abstract,
                "journal": journal_name,
                "year": year,
            }
        except Exception as parse_err:
            print(f"  ⚠ 解析某 PMID 失败: {parse_err}")
    
    return result


def analyze_abstract(abstract_text, title=""):
    """
    分析摘要文本，判断是否含引物序列、标准曲线、Ct 表格。
    返回: (has_primer_seq, has_std_curve, has_ct_table)
    """
    text = (title + " " + abstract_text).lower()
    
    # === has_primer_seq: 检测引物/探针序列特征 ===
    primer_patterns = [
        r"forward\s*primer.*?(?:5[\'\u2019]?\s*[:-]?\s*([atcg]{15,}))",
        r"reverse\s*primer.*?(?:5[\'\u2019]?\s*[:-]?\s*([atcg]{15,}))",
        r"probe.*?(?:5?[\'\u2019]?\s*[:-]?\s*(?:fam|hex|vic|cy[35]|tet|tamra|bhq|mgb).*?([atcg]{10,}))",
        r"(?:forward|reverse|probe|primer|fwd|rev).*?sequence.*?(?:were|was|are|is|:)\s*[\'\u2019]?\s*([atcg]{15,})",
        r"primers?\s*(?:used|were|designed).*?(?:fw|rv|fwd|rev|forward|reverse|probe)",
        r"(?:fw|rv|fwd|rev)\s*(?:primer)?\s*[\'\u2019]?\s*(?:sequence)?.*?[atcg]{15,}",
        r"nucleotide\s*sequence.*?(?:of|for).*?primer",
        r"oligonucleotide\s*(?:primer|probe)",
        r"(?:table|listed|shown).*?(?:primer|probe).*?(?:sequence|oligo)",
        r"primers?\s+and\s+probes?\s+(?:were|are|used|designed|listed)",
    ]
    has_primer_seq = False
    for pat in primer_patterns:
        if re.search(pat, text, re.IGNORECASE):
            has_primer_seq = True
            break
    
    # === has_std_curve: 检测标准曲线参数 ===
    std_curve_patterns = [
        r"standard\s*curve",
        r"slope\s*(?:of|was|=)\s*(-?\d+\.?\d*)",
        r"intercept\s*(?:of|was|=)\s*(\d+\.?\d*)",
        r"r\s*(?:square|squared|2|²)\s*(?:of|was|=|>|≥|:)\s*(0?\.\d+)",
        r"efficiency\s*(?:of|was|=)\s*(\d+\.?\d*)\s*%",
        r"pcr\s*efficiency",
        r"amplification\s*efficiency",
        r"correlation\s*coefficient",
        r"(?:y\s*=\s*-?\d+\.?\d*\s*x\s*\+\s*\d+\.?\d*)",  # 回归方程
        r"limit\s*of\s*detection.*?(?:copies|copy|cp)",
        r"lod.*?(?:copies|copy|cp|reaction)",
        r"detection\s*limit.*?copies",
        r"calibration\s*curve",
        r"(?:log|log10).*?(?:copies|copy).*?(?:ct|ct value|cq)",
        r"(?:copies|copy)\s*(?:number|per).*?(?:ct|ct value|cq)",
    ]
    has_std_curve = False
    for pat in std_curve_patterns:
        if re.search(pat, text, re.IGNORECASE):
            has_std_curve = True
            break
    
    # === has_ct_table: 检测是否含逐点 Ct 数据 ===
    ct_patterns = [
        r"(?:ct|ct\s*value|cq|cq\s*value|threshold\s*cycle).*?(?:table|data|shown|listed|summar)",
        r"mean\s*(?:ct|ct\s*value|cq).*?(?:±|sd|standard\s*deviation)",
        r"(?:table|fig).*?(?:ct|ct\s*value|cq|cq\s*value|threshold\s*cycle)",
    ]
    has_ct_table = False
    for pat in ct_patterns:
        if re.search(pat, text, re.IGNORECASE):
            has_ct_table = True
            break
    
    return has_primer_seq, has_std_curve, has_ct_table


def main():
    print("=" * 70)
    print("qPCR 8病原体 × 3轮 PubMed 文献检索")
    print(f"开始时间: {datetime.now().isoformat()}")
    print("=" * 70)
    
    all_results = []
    seen_pmids = set()
    stats = {}
    
    total_searches = sum(len(p["queries"]) for p in PATHOGENS)
    search_idx = 0
    
    for pathogen in PATHOGENS:
        pname = pathogen["name"]
        print(f"\n{'─' * 70}")
        print(f"🔬 病原体: {pname} ({pathogen['genes']})")
        print(f"{'─' * 70}")
        
        for round_num in sorted(pathogen["queries"].keys()):
            search_idx += 1
            query = pathogen["queries"][round_num]
            print(f"\n  [R{round_num}] 检索 {search_idx}/{total_searches}")
            print(f"  Query: {query[:120]}...")
            
            # esearch
            esearch_results = safe_esearch(query)
            if not esearch_results:
                print(f"  ❌ 检索失败，跳过")
                stats.setdefault(pname, {})[f"R{round_num}"] = 0
                time.sleep(SLEEP_BETWEEN_REQUESTS)
                continue
            
            all_pmids = esearch_results.get("IdList", [])
            total_count = int(esearch_results.get("Count", 0))
            print(f"  📊 PubMed 命中: {total_count}，返回 {len(all_pmids)} 篇")
            
            if not all_pmids:
                stats.setdefault(pname, {})[f"R{round_num}"] = 0
                time.sleep(SLEEP_BETWEEN_REQUESTS)
                continue
            
            # efetch
            time.sleep(SLEEP_BETWEEN_REQUESTS)
            abstracts = safe_efetch(all_pmids)
            
            round_results = 0
            for pmid in all_pmids:
                if pmid not in abstracts:
                    continue
                info = abstracts[pmid]
                
                has_primer, has_curve, has_ct = analyze_abstract(
                    info["abstract"], info["title"]
                )
                
                # 筛选：保留至少含引物序列或标准曲线的文献
                # 但为完整性保留所有记录，用标志位标记
                entry = {
                    "pmid": int(pmid),
                    "title": info["title"],
                    "journal": info["journal"],
                    "year": info["year"],
                    "pathogen": pname,
                    "search_round": round_num,
                    "has_primer_seq": has_primer,
                    "has_std_curve": has_curve,
                    "has_ct_table": has_ct,
                    "notes": f"genes: {pathogen['genes']}",
                }
                
                if pmid not in seen_pmids:
                    seen_pmids.add(pmid)
                    all_results.append(entry)
                    round_results += 1
                else:
                    # 已存在但更新 search_round 标记
                    for existing in all_results:
                        if existing["pmid"] == int(pmid):
                            existing["search_round"] = max(existing["search_round"], round_num)
                            # 合并标志位（OR 逻辑）
                            existing["has_primer_seq"] = existing["has_primer_seq"] or has_primer
                            existing["has_std_curve"] = existing["has_std_curve"] or has_curve
                            existing["has_ct_table"] = existing["has_ct_table"] or has_ct
                            break
            
            stats.setdefault(pname, {})[f"R{round_num}"] = round_results
            print(f"  ✅ 新增 {round_results} 篇（去重后）")
            
            time.sleep(SLEEP_BETWEEN_REQUESTS)
    
    # ========== 输出统计 ==========
    print(f"\n{'=' * 70}")
    print("检索完成统计")
    print(f"{'=' * 70}")
    
    total = len(all_results)
    with_primer = sum(1 for r in all_results if r["has_primer_seq"])
    with_curve = sum(1 for r in all_results if r["has_std_curve"])
    with_ct = sum(1 for r in all_results if r["has_ct_table"])
    with_both = sum(1 for r in all_results if r["has_primer_seq"] and r["has_std_curve"])
    
    print(f"\n总去重文献数: {total}")
    print(f"  含引物序列: {with_primer} ({with_primer/total*100:.1f}%)" if total else "  含引物序列: 0")
    print(f"  含标准曲线: {with_curve} ({with_curve/total*100:.1f}%)" if total else "  含标准曲线: 0")
    print(f"  含Ct数据表: {with_ct} ({with_ct/total*100:.1f}%)" if total else "  含Ct数据表: 0")
    print(f"  同时含引物+标准曲线: {with_both}")
    
    print(f"\n分病原体统计:")
    for pname in sorted(set(r["pathogen"] for r in all_results)):
        precs = [r for r in all_results if r["pathogen"] == pname]
        p_with_primer = sum(1 for r in precs if r["has_primer_seq"])
        p_with_curve = sum(1 for r in precs if r["has_std_curve"])
        p_both = sum(1 for r in precs if r["has_primer_seq"] and r["has_std_curve"])
        print(f"  {pname}: {len(precs)} 篇 (引物:{p_with_primer} 曲线:{p_with_curve} 双全:{p_both})")
    
    # ========== 写入 JSON ==========
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(all_results, f, ensure_ascii=False, indent=2)
    
    print(f"\n✅ 结果已保存到: {OUTPUT_FILE}")
    print(f"文件大小: {os.path.getsize(OUTPUT_FILE):,} bytes")
    print(f"文献条目: {len(all_results)}")
    
    # 同时保存统计信息
    stats_file = OUTPUT_FILE.replace(".json", "_stats.json")
    stats_data = {
        "total": total,
        "with_primer_seq": with_primer,
        "with_std_curve": with_curve,
        "with_ct_table": with_ct,
        "with_both": with_both,
        "by_pathogen": {},
        "timestamp": datetime.now().isoformat(),
    }
    for pname in sorted(set(r["pathogen"] for r in all_results)):
        precs = [r for r in all_results if r["pathogen"] == pname]
        stats_data["by_pathogen"][pname] = {
            "total": len(precs),
            "with_primer": sum(1 for r in precs if r["has_primer_seq"]),
            "with_curve": sum(1 for r in precs if r["has_std_curve"]),
            "with_both": sum(1 for r in precs if r["has_primer_seq"] and r["has_std_curve"]),
        }
    
    with open(stats_file, "w", encoding="utf-8") as f:
        json.dump(stats_data, f, ensure_ascii=False, indent=2)
    
    return all_results


if __name__ == "__main__":
    main()
