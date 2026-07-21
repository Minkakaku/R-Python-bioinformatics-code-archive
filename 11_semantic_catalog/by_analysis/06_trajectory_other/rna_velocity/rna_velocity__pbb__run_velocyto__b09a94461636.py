#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import subprocess
import sys
import logging


script_dir = "/mnt/NexData_file/home/pbb/auto/velocyto/test/module"

# 被调用脚本
velocyto_script = os.path.join(
    script_dir,
    "run_velocyto.py"
)


# 日志
def setup_logger():

    logger = logging.getLogger("velocyto_pipeline")
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s'
    )
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


# 加载 Excel
def load_data(xls_path):

    try:

        ext = os.path.splitext(xls_path)[-1].lower()
        if ext == ".xls":
            df = pd.read_excel(
                xls_path,
                dtype=str,
                engine="xlrd"
            )

        elif ext == ".xlsx":

            df = pd.read_excel(
                xls_path,
                dtype=str,
                engine="openpyxl"
            )

        else:

            logger.error(
                f"不支持的文件格式: {ext}"
            )

            return None

        required_columns = [
            '项目编号',
            '样本分析名',
            '物种中文名',
            '捕获类型'
        ]

        missing_cols = [
            col for col in required_columns
            if col not in df.columns
        ]

        if missing_cols:

            logger.error(
                f"Excel 缺少列: {', '.join(missing_cols)}"
            )
            return None
        return df

    except Exception as e:
        logger.error(
            f"加载 Excel 失败: {str(e)}"
        )
        return None


# 物种映射
species_map = {

    "人": "human",
    "小鼠": "mouse",
    "大鼠": "rat",
    "人+小鼠": "double",
    "鼠+人": "double",
    "巴马香猪": "bama_pig",
    "新西兰兔": "rabbit",
    "兔子": "rabbit",
    "鸡": "chicken"
}


# 获取 rawdata 目录
def get_base_dir(project_id):

    if project_id.startswith("SC26"):
        return "/mnt/NexData_file/nas/rawdata/rawdata-sc-2026"

    elif project_id.startswith("SC25"):
        return "/mnt/NexData_file/nas/rawdata/rawdata-sc-2025"

    else:
        raise ValueError(
            f"未知项目编号: {project_id}"
        )


# 处理单个样本
def process_sample(row):

    try:
        project_id = str(
            row['项目编号']
        ).split("-")[0]

        sample_name = row['样本分析名']
        species_cn = row['物种中文名']
        capture_type = row['捕获类型']
        logger.info(
            f"开始处理样本: "
            f"{sample_name} "
            f"(项目: {project_id}, "
            f"物种: {species_cn}, "
            f"捕获类型: {capture_type})"
        )

        # 物种转换
        species = species_map.get(species_cn)

        if not species:
            logger.error(
                f"不支持的物种: "
                f"{species_cn}"
            )
            return

        # 只处理支持的物种
        if species not in ["human", "mouse"]:
            logger.warning(
                f"velocyto 暂不支持物种: "
                f"{species}"
            )
            return

        # 获取项目目录
        base_dir = get_base_dir(project_id)
        sample_path = os.path.join(
            base_dir,
            project_id,
            sample_name
        )

        # 检查样本目录
        if not os.path.exists(sample_path):
            logger.error(
                f"样本目录不存在: "
                f"{sample_path}"
            )

            return

        # velocyto 输出目录
        velocyto_dir = os.path.join(
            sample_path,
            "velocyto"
        )

        os.makedirs(
            velocyto_dir,
            exist_ok=True
        )

        # loom 文件
        loom_file = os.path.join(
            velocyto_dir,
            f"{sample_name}.loom"
        )

        # 已完成则跳过
        if os.path.exists(loom_file):

            logger.warning(
                f"loom 文件已存在，跳过: "
                f"{sample_name}"
            )

            return

        # log
        log_file = os.path.join(
            velocyto_dir,
            f"{sample_name}.sbatch.log"
        )

        python_env = (
            "/mnt/NexData_file/share/scRNA-rawdata/envs/scRawdata/bin/python"
        )

        # sbatch 命令
        cmd = (
            f"sbatch "
            f"--job-name={sample_name}_velocyto "
            f"--partition=all "
            f"--cpus-per-task=32 "
            f"--mem=256G "
            f"--time=48:00:00 "
            f"--output={log_file} "
            f"--wrap='"
            f"{python_env} "
            f"{velocyto_script} "
            f"{sample_path} "
            f"{species} "
            f"\"{capture_type}\""
            f"'"
        )

        logger.info(f"执行命令: {cmd}")

        subprocess.run(
            cmd,
            shell=True,
            check=True
        )

    except Exception as e:

        logger.error(
            f"处理样本 "
            f"{row.get('样本分析名', '未知')} "
            f"时出错: "
            f"{str(e)}"
        )


# 主函数
def main():

    parser = argparse.ArgumentParser(
        description='velocyto 自动提交脚本'
    )

    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help='输入 Excel 文件'
    )

    args = parser.parse_args()
    global logger
    logger = setup_logger()
    logger.info("启动 velocyto 流程")
    logger.info(
        f"输入文件: {args.input}"
    )

    # 加载数据
    df = load_data(args.input)
    if df is None:
        logger.error("读取 Excel 失败")
        sys.exit(1)
    logger.info(
        f"共读取 "
        f"{len(df)} 个样本"
    )

    # 循环提交
    for _, row in df.iterrows():
        process_sample(row)
    logger.info("全部任务提交完成")

if __name__ == "__main__":

    main()