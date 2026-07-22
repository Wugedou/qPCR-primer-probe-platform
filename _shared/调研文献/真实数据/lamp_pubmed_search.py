#!/usr/bin/env python3
"""
LAMP 8 Pathogen × 3 Round PubMed Literature Search
Output: lamp_search_raw.json
"""
import json
import time
import sys
from datetime import datetime
from Bio import Entrez

# ─── CONFIG ───────────────────────────────────────────────────────
Entrez.email = "zhangxiaoguang@example.com"  # NCBI requires email
MAX_RESULTS = 20      # top N per search round
SLEEP = 0.5           # seconds between API calls (NCBI: ≤3/sec without API key)

OUTPUT_PATH = "/Users/zhangxiaoguang/Documents/16人工智能AI/projects/qPCR_引物探针设计仿真平台/_shared/调研文献/真实数据/lamp_search_raw.json"

# ─── PATHOGENS & SEARCH STRATEGIES ────────────────────────────────
PATHOGENS = {
    "SARS-CoV-2": {
        "name_en": "SARS-CoV-2",
        "name_cn": "新型冠状病毒",
        "taxonomy": "Virus / SARS-CoV-2",
        "round_1": '("SARS-CoV-2" OR "COVID-19") AND LAMP AND (Tt OR "time to positive" OR "standard curve" OR "real-time") AND (copies OR quantification)',
        "round_2": '("SARS-CoV-2" OR "COVID-19") AND LAMP AND (primer OR F3 OR B3 OR FIP OR BIP) AND sequence',
        "round_3": '("SARS-CoV-2" OR "COVID-19") AND LAMP AND (specificity OR cross-reactivity OR LOD)',
    },
    "Influenza_A": {
        "name_en": "Influenza A virus",
        "name_cn": "甲型流感病毒",
        "taxonomy": "Virus / Influenza A",
        "round_1": '("Influenza A" OR "IAV") AND LAMP AND (Tt OR "time to positive" OR "standard curve" OR "real-time") AND (copies OR quantification)',
        "round_2": '("Influenza A" OR "IAV") AND LAMP AND (primer OR F3 OR B3 OR FIP OR BIP) AND sequence',
        "round_3": '("Influenza A" OR "IAV") AND LAMP AND (specificity OR cross-reactivity OR LOD)',
    },
    "Dengue": {
        "name_en": "Dengue virus",
        "name_cn": "登革病毒",
        "taxonomy": "Virus / Dengue",
        "round_1": '(Dengue OR DENV) AND LAMP AND (Tt OR "time to positive" OR "standard curve" OR "real-time") AND (copies OR quantification)',
        "round_2": '(Dengue OR DENV) AND LAMP AND (primer OR F3 OR B3 OR FIP OR BIP) AND sequence',
        "round_3": '(Dengue OR DENV) AND LAMP AND (specificity OR cross-reactivity OR LOD)',
    },
    "HBV": {
        "name_en": "Hepatitis B virus",
        "name_cn": "乙型肝炎病毒",
        "taxonomy": "Virus / HBV",
        "round_1": '(HBV OR "Hepatitis B") AND LAMP AND (Tt OR "time to positive" OR "standard curve" OR "real-time") AND (copies OR quantification)',
        "round_2": '(HBV OR "Hepatitis B") AND LAMP AND (primer OR F3 OR B3 OR FIP OR BIP) AND sequence',
        "round_3": '(HBV OR "Hepatitis B") AND LAMP AND (specificity OR cross-reactivity OR LOD)',
    },
    "MTB_IS6110": {
        "name_en": "Mycobacterium tuberculosis (IS6110)",
        "name_cn": "结核分枝杆菌 (IS6110)",
        "taxonomy": "Bacteria / MTB",
        "round_1": '("Mycobacterium tuberculosis" OR MTB) AND (IS6110 OR IS1081) AND LAMP AND (Tt OR "time to positive" OR "standard curve" OR "real-time") AND (copies OR quantification)',
        "round_2": '("Mycobacterium tuberculosis" OR MTB) AND (IS6110 OR IS1081) AND LAMP AND (primer OR F3 OR B3 OR FIP OR BIP) AND sequence',
        "round_3": '("Mycobacterium tuberculosis" OR MTB) AND (IS6110 OR IS1081) AND LAMP AND (specificity OR cross-reactivity OR LOD)',
    },
    "HIV-1": {
        "name_en": "HIV-1",
        "name_cn": "人类免疫缺陷病毒1型",
        "taxonomy": "Virus / HIV-1",
        "round_1": '("HIV-1" OR HIV) AND LAMP AND (Tt OR "time to positive" OR "standard curve" OR "real-time") AND (copies OR quantification)',
        "round_2": '("HIV-1" OR HIV) AND LAMP AND (primer OR F3 OR B3 OR FIP OR BIP) AND sequence',
        "round_3": '("HIV-1" OR HIV) AND LAMP AND (specificity OR cross-reactivity OR LOD)',
    },
    "ZIKV": {
        "name_en": "Zika virus",
        "name_cn": "寨卡病毒",
        "taxonomy": "Virus / Zika",
        "round_1": '(Zika OR ZIKV) AND LAMP AND (Tt OR "time to positive" OR "standard curve" OR "real-time") AND (copies OR quantification)',
        "round_2": '(Zika OR ZIKV) AND LAMP AND (primer OR F3 OR B3 OR FIP OR BIP) AND sequence',
        "round_3": '(Zika OR ZIKV) AND LAMP AND (specificity OR cross-reactivity OR LOD)',
    },
    "E_coli_O157": {
        "name_en": "E. coli O157:H7",
        "name_cn": "大肠杆菌 O157:H7",
        "taxonomy": "Bacteria / E. coli O157:H7",
        "round_1": '("E. coli O157" OR "Escherichia coli O157") AND LAMP AND (Tt OR "time to positive" OR "standard curve" OR "real-time") AND (copies OR quantification)',
        "round_2": '("E. coli O157" OR "Escherichia coli O157") AND LAMP AND (primer OR F3 OR B3 OR FIP OR BIP) AND sequence',
        "round_3": '("E. coli O157" OR "Escherichia coli O157") AND LAMP AND (specificity OR cross-reactivity OR LOD)',
    },
}

FOUNDATIONAL = {
    "Notomi_2000": {"pmid": "10807277", "search_extra": False},
    "Tomita_2008": {
        "pmid": None,
        "search_query": 'Tomita AND (LAMP OR "loop-mediated isothermal") AND (2008[dp])',
        "search_extra": True,
    },
    "Nagamine_2002": {
        "pmid": None,
        "search_query": 'Nagamine AND (LAMP OR "loop-mediated isothermal") AND (2002[dp]) AND (loop primer OR "loop primers" OR "accelerated")',
        "search_extra": True,
    },
}


def safe_entrez(func, max_retries=3, **kwargs):
    """Wrapper with retry for NCBI rate limits."""
    for attempt in range(max_retries):
        try:
            time.sleep(SLEEP)
            return func(**kwargs)
        except Exception as e:
            if attempt == max_retries - 1:
                print(f"  ⚠ API error after {max_retries} retries: {e}", file=sys.stderr)
                return None
            wait = (attempt + 1) * 5
            print(f"  ↻ Retry {attempt+1}/{max_retries} after {wait}s...", file=sys.stderr)
            time.sleep(wait)


def search_pubmed(query, max_results=MAX_RESULTS):
    """Run a PubMed search and return list of PMIDs."""
    print(f"  🔍 Search: {query[:120]}...", file=sys.stderr)
    handle = safe_entrez(Entrez.esearch, db="pubmed", term=query,
                         retmax=max_results, sort="relevance")
    if handle is None:
        return []
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    count = int(record.get("Count", "0"))
    print(f"     → Found {count} total, fetched {len(ids)} PMIDs", file=sys.stderr)
    return ids


def fetch_details(pmids):
    """Fetch article details for a list of PMIDs."""
    if not pmids:
        return []
    print(f"  📥 Fetching {len(pmids)} records...", file=sys.stderr)
    handle = safe_entrez(Entrez.efetch, db="pubmed", id=",".join(pmids),
                         rettype="xml", retmode="xml")
    if handle is None:
        return []
    records = Entrez.read(handle)
    handle.close()
    return records.get("PubmedArticle", [])


def parse_article(article_xml):
    """Extract structured info from a PubmedArticle XML element."""
    medline = article_xml.get("MedlineCitation", {})
    art = medline.get("Article", {})

    pmid = str(medline.get("PMID", ""))
    
    # Title
    title = art.get("ArticleTitle", "")
    
    # Authors
    authors = []
    author_list = art.get("AuthorList", [])
    if isinstance(author_list, list) and author_list:
        for a in author_list[:10]:  # first 10 authors
            try:
                last = a.get("LastName", "")
                fore = a.get("ForeName", "")
                if last:
                    authors.append(f"{last} {fore}".strip())
            except Exception:
                pass
    
    # Journal
    journal_info = art.get("Journal", {})
    journal = journal_info.get("Title", "")
    isoabbr = journal_info.get("ISOAbbreviation", "")
    
    # Year
    journal_issue = journal_info.get("JournalIssue", {})
    pub_date = journal_issue.get("PubDate", {})
    year = pub_date.get("Year", "")
    if not year:
        medline_date = pub_date.get("MedlineDate", "")
        if medline_date:
            year = medline_date[:4]
    
    # Abstract
    abstract_obj = art.get("Abstract", {})
    abstract_texts = abstract_obj.get("AbstractText", [])
    abstract = ""
    if isinstance(abstract_texts, list):
        parts = []
        for at in abstract_texts:
            label = at.attributes.get("Label", "") if hasattr(at, 'attributes') else ""
            text = str(at)
            if label:
                parts.append(f"{label}: {text}")
            else:
                parts.append(text)
        abstract = " ".join(parts)
    elif isinstance(abstract_texts, str):
        abstract = abstract_texts
    
    # DOI
    doi = ""
    eloc_id = ""
    article_ids = art.get("ELocationID", [])
    if isinstance(article_ids, list):
        for eid in article_ids:
            eid_str = str(eid)
            if eid.attributes.get("EIdType") == "doi":
                doi = eid_str
            elif eid.attributes.get("EIdType") == "pii":
                eloc_id = eid_str
    
    # Publication type
    pub_types = []
    pub_type_list = art.get("PublicationTypeList", [])
    if isinstance(pub_type_list, list):
        for pt in pub_type_list:
            pub_types.append(str(pt))

    # Mesh terms
    mesh_terms = []
    mesh_list = medline.get("MeshHeadingList", [])
    if mesh_list:
        for mh in mesh_list:
            desc = mh.get("DescriptorName", "")
            if desc:
                mesh_terms.append(str(desc))

    # Keywords
    keywords = []
    kw_list = medline.get("KeywordList", [])
    if kw_list:
        for kws in kw_list:
            for kw in kws:
                keywords.append(str(kw))
    
    # ── Article Content Analysis ──────────────────────────────
    full_text = (title + " " + abstract).lower()
    
    # 引物检测 (含6条LAMP引物: F3/B3/FIP/BIP/LF/LB)
    primer_keywords = ["f3", "b3", "fip", "bip", "lf", "lb", "loop primer",
                       "loop f", "loop b", "forward inner", "backward inner",
                       "forward outer", "backward outer", "primers:", "primer set"]
    has_primers = any(kw in full_text for kw in primer_keywords)
    
    # 引物序列检测
    primer_seq_indicators = [
        "F3:", "B3:", "FIP:", "BIP:", "LF:", "LB:",
        "5'-", "5′", "F3 sequence", "B3 sequence",
        "forward outer:", "backward outer:",
        "forward inner:", "backward inner:"
    ]
    has_primer_sequences = any(ind in full_text for ind in primer_seq_indicators)
    
    # Tt / time-to-threshold 数据
    tt_indicators = ["tt", "time to positive", "time-to-positive", "threshold time",
                     "time to threshold", "reaction time", "amplification time",
                     "standard curve", "calibration curve", "copy number"]
    has_tt_data = any(ind in full_text for ind in tt_indicators)
    
    # LOD 数据
    lod_indicators = ["lod", "limit of detection", "detection limit",
                      "sensitivity", "analytical sensitivity", "loq",
                      "limit of quantification", "copies per", "copy/",
                      "cfu/ml", "pfu/ml", "copies/μl", "copies/µl",
                      "minimum detection"]
    has_lod_data = any(ind in full_text for ind in lod_indicators)

    # 引物条数检测
    primer_count = None
    for indicator in ["6 primers", "6-primer", "six primers", "six-primer", "six primer",
                       "5 primers", "5-primer", "five primers", "4 primers", "4-primer", "four primers"]:
        if indicator in full_text:
            primer_count = indicator
            break

    # Real-time / 定量检测
    is_real_time = any(kw in full_text for kw in 
                       ["real-time", "real time", "quantitative", "intercalat",
                        "fluorescence", "fluorogenic", "sybr", "syto", "eva green",
                        "calcein", "turbidimetry", "turbidity", "colorimetric"])
    
    return {
        "pmid": pmid,
        "title": title,
        "authors": authors,
        "journal": journal or isoabbr,
        "year": year,
        "doi": doi,
        "abstract": abstract[:1500] if abstract else "",  # truncate long abstracts
        "publication_types": pub_types,
        "mesh_terms": mesh_terms[:10],
        "keywords": keywords[:10],
        # Analysis flags
        "has_primers": has_primers,
        "has_primer_sequences": has_primer_sequences,
        "primer_count_hint": primer_count,
        "has_tt_data": has_tt_data,
        "has_lod_data": has_lod_data,
        "is_real_time": is_real_time,
    }


def search_one_round(pathogen_key, round_num, query, total_hits_known):
    """Execute one round of search and return parsed results."""
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"📌 {pathogen_key} — Round {round_num}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    
    pmids = search_pubmed(query)
    if not pmids:
        return {"query": query, "total_hits_pubmed": total_hits_known, "results_fetched": 0, "results": []}
    
    articles = fetch_details(pmids)
    results = []
    for art in articles:
        try:
            parsed = parse_article(art)
            results.append(parsed)
        except Exception as e:
            print(f"  ⚠ Parse error for article: {e}", file=sys.stderr)
    
    return {
        "query": query,
        "total_hits_pubmed": total_hits_known,
        "results_fetched": len(results),
        "results": results,
    }


def main():
    print("╔══════════════════════════════════════════════════════════╗", file=sys.stderr)
    print("║  LAMP 8 Pathogen × 3 Round PubMed Literature Search    ║", file=sys.stderr)
    print("╚══════════════════════════════════════════════════════════╝", file=sys.stderr)
    print(f"Started: {datetime.now().isoformat()}", file=sys.stderr)
    print(f"Output: {OUTPUT_PATH}", file=sys.stderr)
    
    output = {
        "metadata": {
            "search_date": datetime.now().strftime("%Y-%m-%d"),
            "search_iso": datetime.now().isoformat(),
            "database": "PubMed (via Bio.Entrez)",
            "max_results_per_round": MAX_RESULTS,
            "total_pathogens": len(PATHOGENS),
            "rounds_per_pathogen": 3,
            "note": "筛选条件：含6条引物序列(至少F3/B3/FIP/BIP) + 含Tt或LOD数据。排除无引物序列/纯定性文献。",
        },
        "pathogens": {},
        "foundational_papers": {},
    }
    
    # ── Phase 1: 8 Pathogens × 3 Rounds ───────────────────────
    for pkey, pinfo in PATHOGENS.items():
        print(f"\n{'🦠'*30}", file=sys.stderr)
        print(f"🦠 PATHOGEN: {pinfo['name_en']} ({pinfo['name_cn']})", file=sys.stderr)
        
        pathogen_results = {
            "name_en": pinfo["name_en"],
            "name_cn": pinfo["name_cn"],
            "taxonomy": pinfo["taxonomy"],
        }
        
        for round_label in ["round_1", "round_2", "round_3"]:
            round_num = int(round_label.split("_")[1])
            query = pinfo[round_label]
            
            # Quick count first
            count_handle = safe_entrez(Entrez.esearch, db="pubmed", term=query, retmax=1)
            total_count = 0
            if count_handle:
                total_count = int(Entrez.read(count_handle)["Count"])
                count_handle.close()
            
            search_data = search_one_round(pkey, round_num, query, total_count)
            search_data["total_hits_pubmed"] = total_count
            
            # Add round summary
            results_with_flags = search_data.get("results", [])
            search_data["summary"] = {
                "total_hits": total_count,
                "fetched": len(results_with_flags),
                "with_primers": sum(1 for r in results_with_flags if r.get("has_primers")),
                "with_primer_seq": sum(1 for r in results_with_flags if r.get("has_primer_sequences")),
                "with_tt_data": sum(1 for r in results_with_flags if r.get("has_tt_data")),
                "with_lod_data": sum(1 for r in results_with_flags if r.get("has_lod_data")),
                "real_time": sum(1 for r in results_with_flags if r.get("is_real_time")),
            }
            
            pathogen_results[f"round_{round_num}"] = search_data
        
        # Merge summaries
        all_results = []
        for r in range(1, 4):
            all_results.extend(pathogen_results[f"round_{r}"].get("results", []))
        
        # Dedup by PMID
        seen = set()
        unique = []
        for r in all_results:
            if r["pmid"] not in seen:
                seen.add(r["pmid"])
                unique.append(r)
        
        pathogen_results["overall_summary"] = {
            "total_unique_pmids": len(unique),
            "with_primers": sum(1 for r in unique if r.get("has_primers")),
            "with_primer_seq": sum(1 for r in unique if r.get("has_primer_sequences")),
            "with_tt_data": sum(1 for r in unique if r.get("has_tt_data")),
            "with_lod_data": sum(1 for r in unique if r.get("has_lod_data")),
            "candidate_papers": sum(1 for r in unique if r.get("has_primers") and r.get("has_tt_data")),
        }
        
        output["pathogens"][pkey] = pathogen_results
    
    # ── Phase 2: Foundational Papers ──────────────────────────
    print(f"\n{'📚'*30}", file=sys.stderr)
    print("📚 Foundational LAMP Papers", file=sys.stderr)
    
    for fkey, finfo in FOUNDATIONAL.items():
        print(f"\n  📄 {fkey}", file=sys.stderr)
        
        if finfo.get("pmid"):
            # Fetch by PMID
            print(f"     PMID: {finfo['pmid']}", file=sys.stderr)
            articles = fetch_details([finfo["pmid"]])
            if articles:
                try:
                    parsed = parse_article(articles[0])
                    output["foundational_papers"][fkey] = {
                        "pmid": finfo["pmid"],
                        "search_method": "direct_pmid",
                        "result": parsed,
                    }
                except Exception as e:
                    output["foundational_papers"][fkey] = {"error": str(e)}
            else:
                output["foundational_papers"][fkey] = {
                    "pmid": finfo["pmid"],
                    "error": "Not found in PubMed"
                }
                print(f"     ⚠ Could not fetch PMID {finfo['pmid']}", file=sys.stderr)
        elif finfo.get("search_extra"):
            pmids = search_pubmed(finfo["search_query"], max_results=5)
            if pmids:
                articles = fetch_details(pmids)
                results = []
                for art in articles:
                    try:
                        results.append(parse_article(art))
                    except Exception as e:
                        results.append({"error": str(e)})
                output["foundational_papers"][fkey] = {
                    "search_query": finfo["search_query"],
                    "search_method": "search_by_author_year",
                    "results": results,
                }
            else:
                output["foundational_papers"][fkey] = {
                    "search_query": finfo["search_query"],
                    "search_method": "search_by_author_year",
                    "results": [],
                    "error": "No results found"
                }
    
    # ── Output ─────────────────────────────────────────────────
    print(f"\n{'💾'*30}", file=sys.stderr)
    print(f"💾 Writing output to {OUTPUT_PATH}...", file=sys.stderr)
    
    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)
    
    # Print summary statistics
    total_articles = sum(
        len(p.get("round_1", {}).get("results", [])) +
        len(p.get("round_2", {}).get("results", [])) +
        len(p.get("round_3", {}).get("results", []))
        for p in output["pathogens"].values()
    )
    
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"✅ COMPLETE", file=sys.stderr)
    print(f"   Pathogens searched: {len(PATHOGENS)}", file=sys.stderr)
    print(f"   Rounds per pathogen: 3", file=sys.stderr)
    print(f"   Total articles fetched: {total_articles}", file=sys.stderr)
    print(f"   Foundational papers: {len(FOUNDATIONAL)}", file=sys.stderr)
    print(f"   Output: {OUTPUT_PATH}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
