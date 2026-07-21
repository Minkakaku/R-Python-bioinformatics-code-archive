import os
import argparse

def generate_sbatch_commands(input_dir, output_dir):
    # 获取输入目录下所有样本
    sample_dirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]

    # 生成每个样本的 sbatch 命令
    sbatch_commands = []
    for sample in sample_dirs:
        # 定义样本相关路径
        in_bam = os.path.join(input_dir, sample, "Analysis", "step2", "STAR", f"{sample}_SortedByCoordinate.bam")
        umi_file = os.path.join(input_dir, sample, "Analysis", "step3", "umi.xls")
        
        # 创建输出目录结构：loom/样本名
        sample_output_dir = os.path.join(output_dir, "loom", sample)
        if not os.path.exists(sample_output_dir):
            os.makedirs(sample_output_dir)
        
        out_bam = os.path.join(sample_output_dir, f"{sample}_addtag.bam")

        # SLURM sbatch 命令
        sbatch_command = f"""
sbatch --job-name={sample} --partition=all --output={os.path.join(output_dir, f"{sample}.log")} --time=24:00:00 --cpus-per-task=32 --mem=256G --wrap="/mnt/NexData_file/share/scRNA-rawdata/software/seeksoultools.1.3.0/seeksoultools utils addtag --inbam {in_bam} --umifile {umi_file} --outbam {out_bam}"
"""
        sbatch_commands.append(sbatch_command)

    return sbatch_commands

def write_sbatch_scripts(sbatch_commands, output_dir):
    # 输出 SLURM 作业脚本到指定的目录
    script_file = os.path.join(output_dir, "batch_submit.sh")
    with open(script_file, 'w') as f:
        f.writelines(sbatch_commands)
    print(f"批量提交命令已写入: {script_file}")

def main(input_dir, output_dir):
    sbatch_commands = generate_sbatch_commands(input_dir, output_dir)
    write_sbatch_scripts(sbatch_commands, output_dir)

if __name__ == "__main__":
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="批量提交 SLURM 作业")
    parser.add_argument('input_dir', help="输入目录路径")
    parser.add_argument('output_dir', help="输出目录路径")
    
    # 解析命令行参数
    args = parser.parse_args()

    # 调用主函数
    main(args.input_dir, args.output_dir)
