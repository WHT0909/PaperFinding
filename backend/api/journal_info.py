import os
import pandas as pd
from pathlib import Path

# ---------------- 数据文件路径 ----------------
DATA_DIR   = Path(__file__).parent.parent / 'database'
CAS_FILE   = DATA_DIR / 'cas_data.csv'

# 确保目录存在
DATA_DIR.mkdir(exist_ok=True)

# ---------------- 加载中科院分区表 ----------------
def _load_cas():
    if CAS_FILE.exists():
        return pd.read_csv(CAS_FILE)
    # 不存在就返回空表
    return pd.DataFrame(columns=['Journal', 'Division'])

_cas_df = _load_cas()

# ---------------- 主函数：仅按期刊名称查询 ----------------
def get_journal_metrics(journal_info: dict) -> dict:
    """
    仅根据期刊名称返回中科院分区（Division）。
    journal_info 必须包含键 "title"。
    返回: {"cas_division": "1区" 或 None}
    """
    title = (journal_info or {}).get("title", "").strip()
    if not title:
        return {"cas_division": None}

    # 预处理期刊名称：转为小写并移除多余空格
    normalized_title = ' '.join(title.split("(")[0].lower().split())
    
    # 预处理CSV中的期刊名称：转为小写并移除多余空格
    normalized_journals = _cas_df['Journal'].apply(lambda x: ' '.join(x.lower().split()))
    
    # 忽略大小写和空格，匹配期刊名称
    match = _cas_df[normalized_journals == normalized_title]
    if match.empty:
        return {"cas_division": None}

    # 取第一条记录的 Division
    return {"cas_division": match.iloc[0]['Division']}