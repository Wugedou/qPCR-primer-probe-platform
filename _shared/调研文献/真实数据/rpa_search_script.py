#!/usr/bin/env python3
"""
RPA 文献检索脚本：8病原体 × 3轮PubMed检索 + 通用方法学检索
输出结构化JSON，格式同qPCR检索结果。

修复：正确处理 Entrez.read() 返回的 PubmedArticleSet 结构。
"""

import json
import time
import re
from datetime import datetime
from Bio import Entrez

# ===== 配置 =====
Entrez.email = "research@example.com"
OUTPUT_PATH = "/Users/zhangxiaoguang/Documents/16人工智能AI/projects/qPCR_引物探针设计仿真平台/_shared/调研文献/真实数据/rpa_search_raw.json"

# ===== 8种病原体 =====
PATHOGENS = {
    "SARS-CoV-2": ('"SARS-CoV-2"[MeSH] OR "SARS-CoV-2" OR "COVID-19"[MeSH] OR COVID-19'),
    "Influenza A": ('"Influenza A virus"[MeSH] OR "Influenza A" OR "H1N1" OR "H3N2" OR "H5N1" OR "H7N9"'),
    "Dengue": ('"Dengue virus"[MeSH] OR Dengue'),
    "HBV": ('"Hepatitis B virus"[MeSH] OR HBV OR "Hepatitis B"'),
    "MTB": ('"Mycobacterium tuberculosis"[MeSH] OR MTB OR "M. tuberculosis" OR "Mycobacterium tuberculosis"'),
    "HIV-1": ('"HIV-1"[MeSH] OR HIV-1 OR "Human immunodeficiency virus 1"'),
    "ZIKV": ('"Zika virus"[MeSH] OR ZIKV OR "Zika virus"'),
    "E. coli O157": ('"Escherichia coli O157"[MeSH] OR "E. coli O157" OR "EHEC"'),
}

ROUNDS = {
    "R1": ' AND (RPA OR "recombinase polymerase amplification") AND (Tt OR "threshold time" OR "standard curve" OR LOD)',
    "R2": ' AND (RPA OR "recombinase polymerase amplification") AND (primer OR probe OR exo OR nfo) AND sequence',
    "R3": ' AND (RPA OR "recombinase polymerase amplification") AND (specificity OR "cross-reactivity" OR "limit of detection")',
}

GENERAL_QUERY = ' AND (RPA OR "recombinase polymerase amplification") AND ("standard curve" OR quantification) AND (exo OR nfo OR fpg)'


def safe_esearch(query, retmax=30, retry=3):
    for attempt in range(retry):
        try:
            time.sleep(0.4)
            handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, sort="relevance")
            result = Entrez.read(handle)
            handle.close()
            return result
        except Exception as e:
            print(f"  esearch尝试 {attempt+1}/{retry} 失败: {e}")
            if attempt < retry - 1:
                time.sleep(3)
    return None


def safe_efetch(id_list, retry=3):
    """返回 PubmedArticle 列表（从 PubmedArticleSet 中提取）"""
    if not id_list:
        return []
    ids_str = ",".join(id_list)
    for attempt in range(retry):
        try:
            time.sleep(0.5)
            handle = Entrez.efetch(db="pubmed", id=ids_str, rettype="xml", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            # records 是 dict，包含 'PubmedArticle' 键，其值为列表
            articles = records.get("PubmedArticle", [])
            if isinstance(articles, list):
                return articles
            else:
                return [articles] if articles else []
        except Exception as e:
            print(f"  efetch尝试 {attempt+1}/{retry} 失败: {e}")
            if attempt < retry - 1:
                time.sleep(5)
    return []


def safe_str(val, default=""):
    """将 Biopython StringElement 等安全转为 str"""
    try:
        s = str(val)
        return s if s else default
    except:
        return default


def extract_article_info(article_record) -> dict:
    """从单个 PubmedArticle 记录提取信息"""
    medline = article_record.get("MedlineCitation", {})
    article = medline.get("Article", {})

    pmid = safe_str(medline.get("PMID", ""))

    # 标题
    title = safe_str(article.get("ArticleTitle", ""))

    # 期刊
    journal_info = article.get("Journal", {})
    journal_raw = journal_info.get("Title", "")
    journal = safe_str(journal_raw)
    if not journal:
        # 尝试 ISOAbbreviation
        journal = safe_str(journal_info.get("ISOAbbreviation", ""))

    # 年份
    pub_date = journal_info.get("JournalIssue", {}).get("PubDate", {})
    year = safe_str(pub_date.get("Year", ""))
    if not year:
        medline_date = safe_str(pub_date.get("MedlineDate", ""))
        year_match = re.search(r"(\d{4})", medline_date)
        if year_match:
            year = year_match.group(1)

    # 作者
    author_list = article.get("AuthorList", [])
    authors = []
    if isinstance(author_list, list):
        for author in author_list[:5]:
            last = safe_str(author.get("LastName", ""))
            initials = safe_str(author.get("Initials", ""))
            if last:
                authors.append(f"{last} {initials}" if initials else last)

    # 摘要
    abstract_obj = article.get("Abstract", {})
    abstract_parts = abstract_obj.get("AbstractText", [])
    abstract = ""
    if isinstance(abstract_parts, list) and abstract_parts:
        try:
            if isinstance(abstract_parts[0], str):
                abstract = " ".join(safe_str(p) for p in abstract_parts)
            else:
                parts = []
                for part in abstract_parts:
                    label = safe_str(getattr(part, 'attributes', {}).get("Label", ""))
                    text = safe_str(part)
                    parts.append(f"{label}: {text}" if label else text)
                abstract = " ".join(parts)
        except:
            abstract = safe_str(abstract_parts[0])

    # DOI
    doi = ""
    elocation = article.get("ELocationID", [])
    if isinstance(elocation, list):
        for e in elocation:
            eid_type = safe_str(getattr(e, 'attributes', {}).get("EIdType", ""))
            if eid_type == "doi":
                doi = safe_str(e)
                break
        if not doi:
            for e in elocation:
                if "doi" in safe_str(e).lower():
                    doi = safe_str(e)
                    break

    # MeSH 关键词
    mesh_terms = []
    mesh_list = medline.get("MeshHeadingList", [])
    if isinstance(mesh_list, list):
        for mesh in mesh_list:
            desc = safe_str(mesh.get("DescriptorName", ""))
            if desc:
                mesh_terms.append(desc)

    # 关键词
    keywords_list = medline.get("KeywordList", [])
    keywords = []
    if isinstance(keywords_list, list):
        for kw_group in keywords_list:
            if isinstance(kw_group, list):
                for kw in kw_group:
                    keywords.append(safe_str(kw))

    return {
        "PMID": pmid,
        "title": title,
        "journal": journal,
        "year": year,
        "authors": authors,
        "first_author": authors[0] if authors else "",
        "abstract_snippet": abstract[:800] if abstract else "",
        "doi": doi,
        "mesh_terms": mesh_terms[:10],
        "keywords": keywords[:5],
    }


def evaluate_inclusion(article: dict) -> dict:
    """评估文章是否符合筛选标准"""
    title_lower = article["title"].lower()
    abstract_lower = article["abstract_snippet"].lower()
    combined = title_lower + " " + abstract_lower

    has_primer = bool(re.search(r"(primer|forward|reverse|F3|B3|FIP|BIP|LF|LB)", combined))
    has_probe = bool(re.search(r"(probe|exo.probe|nfo.probe|molecular.beacon)", combined))
    has_sequence = bool(re.search(r"(sequence|5'[-\s]|nucleotide)", combined))
    has_tt = bool(re.search(r"(tt|threshold.time|amplification.time)", combined))
    has_lod = bool(re.search(r"(lod|limit.of.detection|sensitivity|detection.limit)", combined))
    has_curve = bool(re.search(r"(standard.curve|calibration.curve|linear.range|linearity|r²|r2)", combined))
    has_probit = bool(re.search(r"(probit)", combined))
    is_lfa = bool(re.search(r"(lateral.flow|immunochromatographic|gold.nanoparticle|colloidal.gold|test.strip|dipstick)", combined))
    is_rpa = bool(re.search(r"(rpa|recombinase.polymerase.amplification)", combined))
    has_fluorescent = bool(re.search(r"(exo|nfo|fpg|fluorescent|real.time|fluorescence)", combined))
    has_real_time = bool(re.search(r"(real.time|qrt.pcr|rt.pcr)", combined))

    tags = []
    if has_primer: tags.append("含引物")
    if has_probe: tags.append("含探针")
    if has_sequence: tags.append("提及序列")
    if has_tt: tags.append("含Tt")
    if has_lod: tags.append("含LOD")
    if has_curve: tags.append("含标准曲线")
    if has_probit: tags.append("含probit")
    if is_lfa: tags.append("LFA胶体金")
    if has_fluorescent: tags.append("荧光检测")
    if has_real_time: tags.append("实时检测")

    exclude_reasons = []
    if not is_rpa:
        exclude_reasons.append("非RPA方法")
    if is_lfa and not has_fluorescent:
        exclude_reasons.append("定性胶体金(非实时荧光)")
    if not has_primer and not has_probe:
        exclude_reasons.append("无引物/探针序列")

    if exclude_reasons:
        status = "excluded"
    elif (has_tt or has_lod or has_curve or has_probit) and (has_primer or has_probe):
        status = "included"
    else:
        status = "candidate"

    return {"tags": tags, "status": status, "exclude_reasons": exclude_reasons}


def run_single_search(label: str, pathogen_query: str, round_key: str, round_query: str, retmax=30) -> list:
    """执行单次检索，返回文章列表（已评估）"""
    full_query = f"({pathogen_query}){round_query}"
    print(f"\n{'='*60}")
    print(f"检索: {label} | {round_key}")
    print(f"Query: {full_query[:130]}...")

    search_result = safe_esearch(full_query, retmax=retmax)
    if not search_result:
        print(f"  ❌ esearch失败")
        return {"articles": [], "count": 0, "total_hits": 0}

    id_list = list(search_result.get("IdList", []))
    total_hits = int(search_result.get("Count", 0))
    print(f"  命中: {total_hits} 篇, 取前 {len(id_list)} 篇")

    if not id_list:
        return {"articles": [], "count": 0, "total_hits": total_hits}

    # 分批efetch (每次最多20篇，NCBI建议)
    articles = []
    batch_size = 20
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        pubmed_articles = safe_efetch(batch)
        for rec in pubmed_articles:
            try:
                info = extract_article_info(rec)
                eval_result = evaluate_inclusion(info)
                info["evaluation"] = eval_result
                articles.append(info)
            except Exception as e:
                print(f"  ⚠️ 解析记录失败: {e}")

    print(f"  解析 {len(articles)} 篇")

    return {"articles": articles, "count": len(articles), "total_hits": total_hits}


def print_summary_line(status_text, articles):
    """打印一行 summary"""
    by_status = {"included": [], "candidate": [], "excluded": []}
    for a in articles:
        by_status.get(a["evaluation"]["status"], []).append(a)
    print(f"  {status_text}: incl={len(by_status['included'])}, "
          f"cand={len(by_status['candidate'])}, excl={len(by_status['excluded'])}")


def main():
    print("=" * 70)
    print("RPA 文献检索系统 v2 — 8病原体 × 3轮 + 通用")
    print("=" * 70)

    all_results = {
        "search_metadata": {
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "tool": "Bio.Entrez / PubMed",
            "email": Entrez.email,
            "target": "RPA Tt标准曲线文献检索 — 引物/探针序列 + Tt/LOD/标准曲线数据",
            "pathogens": list(PATHOGENS.keys()),
            "rounds": {
                "R1": "病原体 + RPA + (Tt OR threshold time OR standard curve OR LOD)",
                "R2": "病原体 + RPA + (primer OR probe OR exo OR nfo) + sequence",
                "R3": "病原体 + RPA + (specificity OR cross-reactivity OR limit of detection)",
            },
            "general_query": 'RPA + ("standard curve" OR quantification) + (exo OR nfo OR fpg)',
            "inclusion_criteria": {
                "include": "含引物/探针序列 + 含逐点Tt数据或LOD/probit数据",
                "exclude": "纯定性胶体金、非RPA、无引物序列",
            },
            "retmax_per_round": 30,
        },
        "results": {},
        "summary": {},
    }

    all_pmids = set()
    total_included = 0
    total_candidate = 0
    total_excluded = 0
    total_fetched = 0

    # === 病原体检索 ===
    for pathogen_key, pathogen_query in PATHOGENS.items():
        all_results["results"][pathogen_key] = {}
        print(f"\n{'#' * 70}")
        print(f"### {pathogen_key}")

        for round_key, round_query in ROUNDS.items():
            result = run_single_search(pathogen_key, pathogen_query, round_key, round_query)
            all_results["results"][pathogen_key][round_key] = {
                "query": f"({pathogen_query}){round_query}",
                "total_hits": result["total_hits"],
                "articles_count": result["count"],
                "articles": result["articles"],
            }
            for a in result["articles"]:
                all_pmids.add(a["PMID"])
                st = a["evaluation"]["status"]
                if st == "included":
                    total_included += 1
                elif st == "candidate":
                    total_candidate += 1
                else:
                    total_excluded += 1
            total_fetched += result["count"]
            print_summary_line("摘要评估", result["articles"])

    # === 通用方法学检索 ===
    print(f"\n{'#' * 70}")
    print("### 通用RPA方法学")

    base_query = '(RPA OR "recombinase polymerase amplification")'
    gen_result = run_single_search("通用方法学", base_query, "方法学", GENERAL_QUERY)
    for a in gen_result["articles"]:
        all_pmids.add(a["PMID"])
        st = a["evaluation"]["status"]
        if st == "included":
            total_included += 1
        elif st == "candidate":
            total_candidate += 1
        else:
            total_excluded += 1
    total_fetched += gen_result["count"]

    all_results["results"]["general_methods"] = {
        "query": f"({base_query}){GENERAL_QUERY}",
        "total_hits": gen_result["total_hits"],
        "articles_count": gen_result["count"],
        "articles": gen_result["articles"],
    }
    print_summary_line("通用方法学", gen_result["articles"])

    # === 汇总 ===
    all_results["summary"] = {
        "total_articles_fetched": total_fetched,
        "total_unique_pmids": len(all_pmids),
        "included_count": total_included,
        "candidate_count": total_candidate,
        "excluded_count": total_excluded,
    }

    # === 写入JSON ===
    print(f"\n{'=' * 70}")
    print("写入输出文件...")

    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(all_results, f, ensure_ascii=False, indent=2)

    # === 最终报告 ===
    s = all_results["summary"]
    print(f"输出: {OUTPUT_PATH}")
    print(f"\n{'=' * 70}")
    print("检索汇总")
    print(f"{'=' * 70}")
    print(f"  检索次数: 25 (8×3 + 1)")
    print(f"  获取篇数(含重复): {s['total_articles_fetched']}")
    print(f"  去重唯一PMID:    {s['total_unique_pmids']}")
    print(f"  保留(included):  {s['included_count']}")
    print(f"  候选(candidate): {s['candidate_count']}")
    print(f"  排除(excluded):  {s['excluded_count']}")

    # 列出保留/候选文献
    print(f"\n{'=' * 70}")
    print("保留 & 候选文献列表")
    print(f"{'=' * 70}")

    for pathogen_key in PATHOGENS:
        for round_key in ROUNDS:
            articles = all_results["results"][pathogen_key][round_key]["articles"]
            kept = [a for a in articles
                    if a["evaluation"]["status"] in ("included", "candidate")]
            if kept:
                print(f"\n--- {pathogen_key} | {round_key} ---")
                for a in kept:
                    tags_str = " | ".join(a["evaluation"]["tags"])
                    print(f"  PMID:{a['PMID']} ({a['year']}) [{a['evaluation']['status']}]")
                    print(f"    {a['title'][:100]}")
                    print(f"    {a['journal']} | 第一作者: {a['first_author']}")
                    print(f"    标签: {tags_str}")

    gen_a = all_results["results"]["general_methods"]["articles"]
    gen_kept = [a for a in gen_a if a["evaluation"]["status"] in ("included", "candidate")]
    if gen_kept:
        print(f"\n--- 通用方法学 ---")
        for a in gen_kept:
            tags_str = " | ".join(a["evaluation"]["tags"])
            print(f"  PMID:{a['PMID']} ({a['year']}) [{a['evaluation']['status']}]")
            print(f"    {a['title'][:100]}")
            print(f"    {a['journal']} | 第一作者: {a['first_author']}")

    print(f"\n✅ 检索完成！")

    # 写一个简洁的markdown版报告
    md_path = OUTPUT_PATH.replace(".json", "_summary.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(f"# RPA Tt标准曲线文献检索结果\n\n")
        f.write(f"**检索日期**: {datetime.now().strftime('%Y-%m-%d')}\n\n")
        f.write(f"## 汇总\n\n")
        f.write(f"| 指标 | 数值 |\n")
        f.write(f"|:---|---:|\n")
        f.write(f"| 检索次数 | 25 |\n")
        f.write(f"| 获取篇数 | {s['total_articles_fetched']} |\n")
        f.write(f"| 去重唯一PMID | {s['total_unique_pmids']} |\n")
        f.write(f"| ✅ 保留 | {s['included_count']} |\n")
        f.write(f"| 🔶 候选 | {s['candidate_count']} |\n")
        f.write(f"| ❌ 排除 | {s['excluded_count']} |\n\n")
        f.write(f"## 保留/候选文献\n\n")
        for pathogen_key in PATHOGENS:
            for round_key in ROUNDS:
                articles = all_results["results"][pathogen_key][round_key]["articles"]
                kept = [a for a in articles
                        if a["evaluation"]["status"] in ("included", "candidate")]
                if kept:
                    f.write(f"### {pathogen_key} | {round_key}\n\n")
                    for a in kept:
                        f.write(f"- **PMID:{a['PMID']}** ({a['year']}) [{a['evaluation']['status']}] — {a['title'][:120]}\n")
                        f.write(f"  - {a['journal']} | {a['first_author']}\n")
                        f.write(f"  - 标签: {' | '.join(a['evaluation']['tags'])}\n")
                    f.write("\n")
        f.write(f"\n完整数据见: `rpa_search_raw.json`\n")

    print(f"Markdown报告: {md_path}")
    return all_results


if __name__ == "__main__":
    main()
