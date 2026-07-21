#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
零宽空格清理工具
用于清理文件中的<200b>（零宽空格）字符
"""

import os
import sys
import re
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Union
import warnings
warnings.filterwarnings('ignore')

class ZeroWidthSpaceCleaner:
    """零宽空格清理器"""
    
    # 零宽空格及相关控制字符的定义
    ZERO_WIDTH_CHARS = {
        '\u200b': 'ZERO WIDTH SPACE (U+200B)',
        '\u200c': 'ZERO WIDTH NON-JOINER (U+200C)',
        '\u200d': 'ZERO WIDTH JOINER (U+200D)',
        '\u200e': 'LEFT-TO-RIGHT MARK (U+200E)',
        '\u200f': 'RIGHT-TO-LEFT MARK (U+200F)',
        '\ufeff': 'BYTE ORDER MARK (U+FEFF)',
        '\u202a': 'LEFT-TO-RIGHT EMBEDDING (U+202A)',
        '\u202b': 'RIGHT-TO-LEFT EMBEDDING (U+202B)',
        '\u202c': 'POP DIRECTIONAL FORMATTING (U+202C)',
        '\u202d': 'LEFT-TO-RIGHT OVERRIDE (U+202D)',
        '\u202e': 'RIGHT-TO-LEFT OVERRIDE (U+202E)',
    }
    
    def __init__(self, verbose: bool = True):
        """
        初始化清理器
        
        参数:
        verbose: 是否显示详细信息
        """
        self.verbose = verbose
        self.stats = {
            'files_processed': 0,
            'total_zwsp_found': 0,
            'files_with_zwsp': 0,
            'cleaned_files': 0
        }
    
    def print_status(self, message: str, indent: int = 0):
        """打印状态信息"""
        if self.verbose:
            indent_str = "  " * indent
            print(f"{indent_str}{message}")
    
    def detect_zwsp_in_string(self, text: str, max_samples: int = 3) -> Dict:
        """
        检测字符串中的零宽空格字符
        
        参数:
        text: 要检测的字符串
        max_samples: 最大样本数
        
        返回:
        检测结果字典
        """
        if not isinstance(text, str):
            return {
                'has_zwsp': False,
                'chars_found': [],
                'samples': []
            }
        
        result = {
            'has_zwsp': False,
            'chars_found': [],
            'samples': []
        }
        
        for char, description in self.ZERO_WIDTH_CHARS.items():
            if char in text:
                result['has_zwsp'] = True
                result['chars_found'].append(description)
                
                # 获取样本
                positions = [m.start() for m in re.finditer(re.escape(char), text)]
                for pos in positions[:max_samples]:
                    start = max(0, pos - 5)
                    end = min(len(text), pos + 5)
                    context = text[start:end]
                    
                    # 用可见标记替换零宽空格
                    visible_context = context.replace(char, '[ZWSP]')
                    result['samples'].append(f"...{visible_context}...")
        
        return result
    
    def remove_zwsp_from_string(self, text: str) -> str:
        """
        从字符串中移除所有零宽空格字符
        
        参数:
        text: 输入字符串
        
        返回:
        清理后的字符串
        """
        if not isinstance(text, str):
            return text
        
        # 移除所有零宽字符
        for char in self.ZERO_WIDTH_CHARS.keys():
            text = text.replace(char, '')
        
        return text
    
    def analyze_file_zwsp(self, filepath: str, sample_size: int = 5) -> Dict:
        """
        分析文件中的零宽空格字符
        
        参数:
        filepath: 文件路径
        sample_size: 每个列检查的样本数
        
        返回:
        分析结果字典
        """
        self.print_status(f"分析文件: {filepath}")
        
        result = {
            'filepath': filepath,
            'has_zwsp': False,
            'columns_with_zwsp': [],
            'total_zwsp_count': 0,
            'zwsp_details': {}
        }
        
        try:
            # 尝试不同的分隔符读取文件
            try:
                df = pd.read_csv(filepath, sep=',')
                separator = ','
            except:
                try:
                    df = pd.read_csv(filepath, sep='\t')
                    separator = '\t'
                except Exception as e:
                    self.print_status(f"  无法读取文件: {e}", indent=1)
                    return result
            
            # 1. 检查列名
            self.print_status("  检查列名...", indent=1)
            for col in df.columns:
                if isinstance(col, str):
                    detection = self.detect_zwsp_in_string(col)
                    if detection['has_zwsp']:
                        result['columns_with_zwsp'].append(col)
                        result['total_zwsp_count'] += len(detection['chars_found'])
                        result['has_zwsp'] = True
            
            # 2. 检查字符串列中的数据
            str_columns = df.select_dtypes(include=['object']).columns.tolist()
            
            for col in str_columns:
                col_has_zwsp = False
                zwsp_samples = []
                
                # 获取样本数据
                for val in df[col].dropna().unique()[:sample_size]:
                    if isinstance(val, str) and val.strip():
                        detection = self.detect_zwsp_in_string(val)
                        if detection['has_zwsp']:
                            col_has_zwsp = True
                            zwsp_samples.append({
                                'value': val,
                                'chars_found': detection['chars_found'],
                                'samples': detection['samples']
                            })
                
                if col_has_zwsp:
                    result['zwsp_details'][col] = {
                        'has_zwsp': True,
                        'sample_count': len(zwsp_samples),
                        'samples': zwsp_samples
                    }
                    result['total_zwsp_count'] += 1
                    result['has_zwsp'] = True
        
        except Exception as e:
            self.print_status(f"  分析文件时出错: {e}", indent=1)
        
        return result
    
    def clean_file(self, input_path: str, output_path: Optional[str] = None, 
                   overwrite: bool = False) -> Tuple[bool, str]:
        """
        清理文件中的零宽空格字符
        
        参数:
        input_path: 输入文件路径
        output_path: 输出文件路径（None则覆盖原文件）
        overwrite: 是否覆盖已存在的输出文件
        
        返回:
        (成功标志, 消息)
        """
        self.stats['files_processed'] += 1
        
        # 确定输出路径
        if output_path is None:
            output_path = input_path
        
        # 检查输出文件是否已存在
        if os.path.exists(output_path) and not overwrite and input_path != output_path:
            return False, f"输出文件已存在: {output_path}"
        
        # 分析文件
        analysis = self.analyze_file_zwsp(input_path)
        
        if not analysis['has_zwsp']:
            self.print_status(f"文件 {input_path} 中没有发现零宽空格字符")
            return True, "没有发现零宽空格字符"
        
        self.stats['files_with_zwsp'] += 1
        
        try:
            # 尝试不同的分隔符读取文件
            try:
                df = pd.read_csv(input_path, sep=',')
                separator = ','
            except:
                try:
                    df = pd.read_csv(input_path, sep='\t')
                    separator = '\t'
                except Exception as e:
                    return False, f"无法读取文件: {e}"
            
            original_shape = df.shape
            zwsp_count = 0
            
            # 1. 清理列名
            original_columns = df.columns.tolist()
            cleaned_columns = []
            
            for col in original_columns:
                if isinstance(col, str):
                    cleaned_col = self.remove_zwsp_from_string(col)
                    cleaned_columns.append(cleaned_col)
                    if cleaned_col != col:
                        zwsp_count += 1
                        self.print_status(f"  清理列名: '{col}' -> '{cleaned_col}'", indent=1)
                else:
                    cleaned_columns.append(col)
            
            df.columns = cleaned_columns
            
            # 2. 清理字符串列中的数据
            str_columns = df.select_dtypes(include=['object']).columns.tolist()
            
            for col in str_columns:
                # 应用清理函数
                df[col] = df[col].apply(self.remove_zwsp_from_string)
                zwsp_count += 1
                self.print_status(f"  清理列: {col}", indent=1)
            
            # 保存清理后的文件
            df.to_csv(output_path, sep=separator, index=False, encoding='utf-8')
            
            self.stats['cleaned_files'] += 1
            
            # 打印摘要
            self.print_status(f"✓ 已清理文件: {input_path}", indent=1)
            self.print_status(f"  原始形状: {original_shape}", indent=2)
            self.print_status(f"  清理后形状: {df.shape}", indent=2)
            self.print_status(f"  发现并清理了 {zwsp_count} 处零宽空格", indent=2)
            
            if analysis['columns_with_zwsp']:
                self.print_status(f"  清理的列名: {analysis['columns_with_zwsp']}", indent=2)
            
            return True, f"成功清理文件，处理了 {zwsp_count} 处零宽空格"
        
        except Exception as e:
            return False, f"清理文件时出错: {e}"
    
    def clean_directory(self, input_dir: str, output_dir: Optional[str] = None, 
                       file_extensions: List[str] = None, recursive: bool = False,
                       overwrite: bool = False) -> Dict:
        """
        清理目录中的所有文件
        
        参数:
        input_dir: 输入目录
        output_dir: 输出目录（None则覆盖原文件）
        file_extensions: 要处理的文件扩展名列表
        recursive: 是否递归处理子目录
        overwrite: 是否覆盖已存在的输出文件
        
        返回:
        处理结果字典
        """
        if file_extensions is None:
            file_extensions = ['.csv', '.tsv', '.txt']
        
        results = {
            'successful': [],
            'failed': [],
            'skipped': []
        }
        
        # 遍历目录
        for root, dirs, files in os.walk(input_dir):
            for file in files:
                # 检查文件扩展名
                if any(file.lower().endswith(ext) for ext in file_extensions):
                    input_path = os.path.join(root, file)
                    
                    # 确定输出路径
                    if output_dir is not None:
                        # 保持目录结构
                        rel_path = os.path.relpath(input_path, input_dir)
                        output_path = os.path.join(output_dir, rel_path)
                        
                        # 创建输出目录
                        os.makedirs(os.path.dirname(output_path), exist_ok=True)
                    else:
                        output_path = input_path
                    
                    # 清理文件
                    success, message = self.clean_file(input_path, output_path, overwrite)
                    
                    if success:
                        results['successful'].append({
                            'input': input_path,
                            'output': output_path,
                            'message': message
                        })
                    else:
                        results['failed'].append({
                            'input': input_path,
                            'output': output_path,
                            'message': message
                        })
            
            # 如果不递归，只处理顶层目录
            if not recursive:
                break
        
        return results
    
    def print_statistics(self):
        """打印统计信息"""
        print("\n" + "="*60)
        print("清理统计")
        print("="*60)
        print(f"处理的文件总数: {self.stats['files_processed']}")
        print(f"包含零宽空格的文件数: {self.stats['files_with_zwsp']}")
        print(f"成功清理的文件数: {self.stats['cleaned_files']}")
        print("="*60)

def main():
    """主函数 - 使用子命令模式"""
    parser = argparse.ArgumentParser(
        description='零宽空格清理工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 清理单个文件
  python remove_zwsp.py single CNs_10_all_cell_features.csv
  
  # 清理单个文件并保存到新文件
  python remove_zwsp.py single input.csv -o cleaned.csv
  
  # 检测文件中的零宽空格
  python remove_zwsp.py detect CNs_10_all_cell_features.csv
  
  # 清理目录中的所有文件
  python remove_zwsp.py directory /path/to/data -o /path/to/cleaned_data
"""
    )
    
    subparsers = parser.add_subparsers(dest='mode', help='运行模式')
    
    # single 模式
    single_parser = subparsers.add_parser('single', help='清理单个文件')
    single_parser.add_argument('input', help='输入文件路径')
    single_parser.add_argument('-o', '--output', help='输出文件路径（默认覆盖原文件）')
    single_parser.add_argument('--overwrite', action='store_true', help='覆盖已存在的输出文件')
    single_parser.add_argument('-q', '--quiet', action='store_true', help='安静模式，不显示详细信息')
    
    # directory 模式
    dir_parser = subparsers.add_parser('directory', help='清理目录中的文件')
    dir_parser.add_argument('input', help='输入目录路径')
    dir_parser.add_argument('-o', '--output', help='输出目录路径（默认覆盖原文件）')
    dir_parser.add_argument('-e', '--extensions', nargs='+', default=['.csv', '.tsv', '.txt'],
                           help='要处理的文件扩展名（默认: .csv .tsv .txt）')
    dir_parser.add_argument('-r', '--recursive', action='store_true', help='递归处理子目录')
    dir_parser.add_argument('--overwrite', action='store_true', help='覆盖已存在的输出文件')
    dir_parser.add_argument('-q', '--quiet', action='store_true', help='安静模式，不显示详细信息')
    
    # detect 模式
    detect_parser = subparsers.add_parser('detect', help='检测文件中的零宽空格字符')
    detect_parser.add_argument('input', help='输入文件路径')
    detect_parser.add_argument('-s', '--sample-size', type=int, default=5, 
                              help='每个列检查的样本数（默认: 5）')
    
    args = parser.parse_args()
    
    # 如果没有提供任何参数，显示帮助
    if not args.mode:
        parser.print_help()
        sys.exit(1)
    
    # 根据模式执行相应的功能
    if args.mode == 'single':
        cleaner = ZeroWidthSpaceCleaner(verbose=not args.quiet)
        
        print(f"清理文件: {args.input}")
        print("-" * 40)
        
        success, message = cleaner.clean_file(args.input, args.output, args.overwrite)
        
        if success:
            print(f"\n✓ {message}")
        else:
            print(f"\n✗ {message}")
            sys.exit(1)
        
        cleaner.print_statistics()
    
    elif args.mode == 'directory':
        cleaner = ZeroWidthSpaceCleaner(verbose=not args.quiet)
        
        print(f"清理目录: {args.input}")
        if args.output:
            print(f"输出目录: {args.output}")
        print(f"文件扩展名: {args.extensions}")
        print(f"递归处理: {args.recursive}")
        print("-" * 40)
        
        results = cleaner.clean_directory(
            args.input, args.output, args.extensions, args.recursive, args.overwrite
        )
        
        # 打印结果摘要
        print(f"\n处理完成:")
        print(f"  成功: {len(results['successful'])} 个文件")
        print(f"  失败: {len(results['failed'])} 个文件")
        print(f"  跳过: {len(results['skipped'])} 个文件")
        
        if results['failed']:
            print("\n失败的文件:")
            for fail in results['failed']:
                print(f"  {fail['input']}: {fail['message']}")
        
        cleaner.print_statistics()
    
    elif args.mode == 'detect':
        cleaner = ZeroWidthSpaceCleaner(verbose=True)
        
        print(f"检测文件: {args.input}")
        print("-" * 40)
        
        analysis = cleaner.analyze_file_zwsp(args.input, args.sample_size)
        
        if analysis['has_zwsp']:
            print(f"\n✗ 发现零宽空格字符!")
            print(f"  总计: {analysis['total_zwsp_count']} 处")
            
            if analysis['columns_with_zwsp']:
                print(f"\n  包含零宽空格的列名:")
                for col in analysis['columns_with_zwsp']:
                    print(f"    - {col}")
            
            if analysis['zwsp_details']:
                print(f"\n  包含零宽空格的数据列:")
                for col, details in analysis['zwsp_details'].items():
                    print(f"    - {col}: {details['sample_count']} 个样本包含")
                    for i, sample in enumerate(details['samples'][:2], 1):  # 只显示前2个样本
                        print(f"        样本 {i}: {sample['value'][:50]}...")
                        for char in sample['chars_found']:
                            print(f"          包含: {char}")
        else:
            print(f"\n✓ 未发现零宽空格字符")

if __name__ == "__main__":
    main()
