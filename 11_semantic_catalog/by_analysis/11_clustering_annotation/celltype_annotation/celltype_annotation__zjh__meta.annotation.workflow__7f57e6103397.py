import argparse
import subprocess
import os
import sys
import multiprocessing # 新增导入
import glob # 新增导入

def run_command(cmd_list, working_dir, error_message_prefix="", log_file_path=None):
    """
    执行一个shell命令并检查其输出，可选择记录到日志文件。
    """
    print(f"Executing command: {' '.join(cmd_list)} in {working_dir}")
    if log_file_path:
        # Ensure log directory exists just before opening the file
        os.makedirs(os.path.dirname(log_file_path), exist_ok=True)
        print(f"Logging to: {log_file_path}")
        with open(log_file_path, 'a', encoding='utf-8') as log_f: # Clear/initialize log for this specific command run
            log_f.write(f"--- EXECUTING COMMAND ---\nCommand: {' '.join(cmd_list)}\nWorking directory: {working_dir}\n\n")

    stdout, stderr = None, None # <--- 在这里初始化
    try:
        if log_file_path:
            with open(log_file_path, 'a', encoding='utf-8') as log_f:
                process = subprocess.Popen(cmd_list, cwd=working_dir, stdout=log_f, stderr=log_f, text=True)
                process.wait() # Wait for the command to complete
        else:
            process = subprocess.Popen(cmd_list, cwd=working_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate() # 如果 log_file_path 为 None，则在这里赋值

        if process.returncode != 0:
            error_detail = f"{error_message_prefix}Error executing command: {' '.join(cmd_list)}\nReturn code: {process.returncode}\n"
            if log_file_path:
                error_detail += f"Check log for details: {log_file_path}\n"
                # Error details are already in the log file from Popen redirection. Add a failure marker.
                with open(log_file_path, 'a', encoding='utf-8') as log_f:
                    log_f.write(f"\n--- COMMAND FAILED ---\nReturn Code: {process.returncode}\n")
            else: # Only print stdout/stderr if not logged to file already
                # 检查stdout是否存在且不为空 (现在 stdout 和 stderr 总是有初始值 None 或被赋值)
                if stdout is not None and stdout.strip(): error_detail += f"Stdout: {stdout}\n"
                # 检查stderr是否存在且不为空
                if stderr is not None and stderr.strip(): 
                    error_detail += f"Stderr: {stderr}\n"
            print(error_detail)
            sys.exit(1)
        else:
            success_message = f"Successfully executed: {' '.join(cmd_list)}"
            print(success_message)
            if log_file_path:
                with open(log_file_path, 'a', encoding='utf-8') as log_f:
                    log_f.write(f"\n--- COMMAND SUCCEEDED ---\n")
            else: # Only print if not logged
                if stdout: print(f"Stdout: {stdout}") # stdout 和 stderr 已被初始化或赋值
                if stderr: print(f"Stderr: {stderr}")
        
        return (None, None) if log_file_path else (stdout, stderr)

    except FileNotFoundError:
        error_msg = f"{error_message_prefix}Error: Command or script not found: {cmd_list[0]}"
        print(error_msg)
        if log_file_path:
            with open(log_file_path, 'a', encoding='utf-8') as log_f:
                log_f.write(f"\n--- FILE NOT FOUND ERROR ---\n{error_msg}\n")
        sys.exit(1)
    except Exception as e:
        error_msg = f"{error_message_prefix}An unexpected error occurred: {e}"
        print(error_msg)
        if log_file_path:
            with open(log_file_path, 'a', encoding='utf-8') as log_f:
                log_f.write(f"\n--- UNEXPECTED ERROR ---\n{error_msg}\n")
        sys.exit(1)

def setup_directories(base_path, subdirs):
    """
    创建目录（如果它们不存在）。
    """
    for subdir in subdirs:
        path = os.path.join(base_path, subdir)
        os.makedirs(path, exist_ok=True)
        print(f"Ensured directory exists: {path}")
    return [os.path.join(base_path, subdir) for subdir in subdirs]

def step1_mode(args, mode):
    """
    执行步骤1：数据预处理和注释（针对特定模式 Pos/Neg）。
    """
    print(f"\n--- Starting Step 1 for {mode} mode ---")
    
    project_data_dir = os.path.join(args.base_project_dir, args.project_name, "data")
    project_analysis_dir = os.path.join(args.base_project_dir, args.project_name, "analysis")
    logs_dir = os.path.join(project_analysis_dir, "logs", mode) # Mode-specific log directory
    os.makedirs(logs_dir, exist_ok=True)
    print(f"Ensured log directory for Step 1 {mode}: {logs_dir}")
    
    mode_data_dir = os.path.join(project_data_dir, mode)
    mode_analysis_dir = os.path.join(project_analysis_dir, mode)
    os.makedirs(mode_analysis_dir, exist_ok=True)

    workflow_script_dir = os.path.join(args.base_project_dir, "workflow_v2", "script")
    workflow_db_dir = os.path.join(args.base_project_dir, "workflow_v2", "database")

    # 1. Preprocessing (Sequential)
    r_preprocessing_script = os.path.join(workflow_script_dir, "preprocessing.no_rt_filter_noAdj.R")
    cmd_preprocess = [
        "Rscript", r_preprocessing_script,
        mode_data_dir, 
        "positive" if mode == "Pos" else "negative"
    ]
    log_file_preprocess = os.path.join(logs_dir, "step1_preprocessing.log")
    run_command(cmd_preprocess, mode_analysis_dir, f"Step 1 ({mode}) Preprocessing: ", log_file_path=log_file_preprocess)

    camera_results_file = os.path.join(mode_analysis_dir, "result_CAMERA.tsv")

    print(f"\n--- Starting Parallel Annotations for {mode} mode ---")
    processes = []
    parallel_commands_info = [] 
    
    # Define commands for parallel execution with log keys
    # 2.1 In-house Annotation
    r_anno_inhouse_script = os.path.join(workflow_script_dir, "anno_inhouse.R")
    db_inhouse = os.path.join(workflow_db_dir, f"in_house_{mode.lower()}.mgf")
    output_inhouse_table = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_inhouse.tsv")
    output_inhouse_table2 = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_inhouse_scores.tsv")
    cmd_anno_inhouse = [
        "Rscript", r_anno_inhouse_script,
        camera_results_file, mode_data_dir, db_inhouse,
        output_inhouse_table, output_inhouse_table2
    ]
    parallel_commands_info.append({"cmd": cmd_anno_inhouse, "name": f"Step 1 ({mode}) In-house Annotation", "log_key": "inhouse"})

    # 2.2 HMDB Annotation
    r_anno_hmdb_script = os.path.join(workflow_script_dir, "anno_HMDB.R")
    db_hmdb1 = os.path.join(workflow_db_dir, "CompDb.Hsapiens.HMDB.5.0.sqlite")
    db_hmdb2 = os.path.join(workflow_db_dir, "CompDb.Hsapiens.HMDB.5.0.2.sqlite")
    ion_file = os.path.join(workflow_db_dir, f"{mode.lower()}_ion.txt")
    ion_type_hmdb = "1" if mode == "Pos" else "0"
    output_hmdb_table = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_HMDB.tsv")
    output_hmdb_table2 = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_HMDB_scores.tsv")
    cmd_anno_hmdb = [
        "Rscript", r_anno_hmdb_script,
        camera_results_file, mode_data_dir, db_hmdb1, db_hmdb2, ion_file, ion_type_hmdb,
        output_hmdb_table, output_hmdb_table2
    ]
    parallel_commands_info.append({"cmd": cmd_anno_hmdb, "name": f"Step 1 ({mode}) HMDB Annotation", "log_key": "hmdb"})

    # 2.3 Lipid Annotation
    r_anno_lipid_script = os.path.join(workflow_script_dir, "anno_lipid.R") 
    db_lipid = os.path.join(workflow_db_dir, f"MSDIAL-TandemMassSpectralAtlas-VS69-{mode}.msp")    
    output_lipid_table = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_lipid.tsv")
    output_lipid_table2 = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_lipid_scores.tsv")
    cmd_anno_lipid = [
        "Rscript", r_anno_lipid_script,
        camera_results_file, mode_data_dir, db_lipid,
        output_lipid_table, output_lipid_table2
    ]
    parallel_commands_info.append({"cmd": cmd_anno_lipid, "name": f"Step 1 ({mode}) Lipid Annotation", "log_key": "lipid"})

    # 2.4 GNPS Annotation
    r_anno_gnps_script = os.path.join(workflow_script_dir, "anno_GNPS.R")  
    db_gnps = os.path.join(workflow_db_dir, f"ALL_GNPS-{mode}.msp")     
    output_gnps_table = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_GNPS.tsv")
    output_gnps_table2 = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_GNPS_scores.tsv")
    cmd_anno_gnps = [
        "Rscript", r_anno_gnps_script,
        camera_results_file, mode_data_dir, db_gnps,
        output_gnps_table, output_gnps_table2
    ]
    parallel_commands_info.append({"cmd": cmd_anno_gnps, "name": f"Step 1 ({mode}) GNPS Annotation", "log_key": "gnps"})

    # 2.5 MoNA Annotation
    r_anno_mona_script = os.path.join(workflow_script_dir, "anno_MoNA.R")  
    db_mona = os.path.join(workflow_db_dir, f"MoNA-export-LC-MS-MS_Spectra-{mode}.msp")     
    output_mona_table = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_MoNA.tsv")
    output_mona_table2 = os.path.join(mode_analysis_dir, "result_CAMERA_annotated_MoNA_scores.tsv")
    cmd_anno_mona = [
        "Rscript", r_anno_mona_script,
        camera_results_file, mode_data_dir, db_mona,
        output_mona_table, output_mona_table2
    ]
    parallel_commands_info.append({"cmd": cmd_anno_mona, "name": f"Step 1 ({mode}) MoNA Annotation", "log_key": "mona"})

    for item in parallel_commands_info:
        cmd = item["cmd"]
        name = item["name"]
        log_key = item["log_key"]
        log_file_name = f"step1_anno_{log_key}.log"
        log_file_path = os.path.join(logs_dir, log_file_name)
        
        print(f"Starting parallel task: {name} - {' '.join(cmd)}. Logging to: {log_file_path}")
        try:
            # Initialize/clear log for this specific command run
            with open(log_file_path, 'w', encoding='utf-8') as log_f: # Use 'w' to overwrite for a new run
                 log_f.write(f"--- EXECUTING COMMAND ---\nCommand: {' '.join(cmd)}\nWorking directory: {mode_analysis_dir}\n\n")
            
            # Append to the log file during Popen
            with open(log_file_path, 'a', encoding='utf-8') as log_f:
                process = subprocess.Popen(cmd, cwd=mode_analysis_dir, stdout=log_f, stderr=log_f, text=True)
            processes.append({"process": process, "name": name, "cmd_str": ' '.join(cmd), "log_file": log_file_path})
        except FileNotFoundError:
            error_msg = f"{name} Error: Command or script not found: {cmd[0]}"
            print(error_msg)
            with open(log_file_path, 'a', encoding='utf-8') as log_f: # Append error to its log
                log_f.write(f"\n--- FILE NOT FOUND ERROR ---\n{error_msg}\n")
            sys.exit(1)
        except Exception as e:
            error_msg = f"{name} An unexpected error occurred during Popen: {e}"
            print(error_msg)
            with open(log_file_path, 'a', encoding='utf-8') as log_f: # Append error to its log
                log_f.write(f"\n--- UNEXPECTED ERROR DURING POPEN ---\n{error_msg}\n")
            sys.exit(1)
            
    for p_info in processes:
        process = p_info["process"]
        name = p_info["name"]
        cmd_str = p_info["cmd_str"]
        log_file = p_info["log_file"]
        
        process.wait() 
        
        if process.returncode != 0:
            error_detail = f"{name} Error executing command: {cmd_str}\nReturn code: {process.returncode}\nCheck log for details: {log_file}"
            print(error_detail)
            with open(log_file, 'a', encoding='utf-8') as log_f:
                log_f.write(f"\n--- COMMAND FAILED ---\nReturn Code: {process.returncode}\n")
            sys.exit(1)
        else:
            print(f"Successfully completed task: {name}. Log: {log_file}")
            with open(log_file, 'a', encoding='utf-8') as log_f:
                log_f.write(f"\n--- COMMAND SUCCEEDED ---\n")
    
    print(f"--- All parallel annotations for {mode} mode finished ---")

    r_anno_inhouse_mod_script = os.path.join(workflow_script_dir, "anno_inhouse_mod.R")
    cmd_anno_inhouse_mod = ["Rscript", r_anno_inhouse_mod_script]
    log_file_inhouse_mod = os.path.join(logs_dir, "step1_anno_inhouse_mod.log")
    run_command(cmd_anno_inhouse_mod, mode_analysis_dir, f"Step 1 ({mode}) In-house Mod Annotation: ", log_file_path=log_file_inhouse_mod)

    py_merge_step1_script = os.path.join(workflow_script_dir, "merge_multiple_annotations_step1.py")
    cmd_merge_step1 = ["python", py_merge_step1_script]
    log_file_merge1 = os.path.join(logs_dir, "step1_merge_annotations_step1.log")
    run_command(cmd_merge_step1, mode_analysis_dir, f"Step 1 ({mode}) Merge Annotations Step 1: ", log_file_path=log_file_merge1)
    
    print(f"--- Finished Step 1 for {mode} mode ---")

def step2_mode(args, mode, num_samples):
    """
    执行步骤2：ClassyFire注释和过滤（针对特定模式 Pos/Neg）。
    """
    print(f"\n--- Starting Step 2 for {mode} mode ---")
    project_analysis_dir = os.path.join(args.base_project_dir, args.project_name, "analysis")
    logs_dir = os.path.join(project_analysis_dir, "logs", mode) # Mode-specific log directory
    os.makedirs(logs_dir, exist_ok=True) # Ensure it exists, might be created by step1
    print(f"Ensured log directory for Step 2 {mode}: {logs_dir}")

    mode_analysis_dir = os.path.join(project_analysis_dir, mode)
    workflow_script_dir = os.path.join(args.base_project_dir, "workflow_v2", "script")

    # 在函数前面统一定义日志文件路径
    log_file_cf_add = os.path.join(logs_dir, "step2_classyfire_add_batch.log")
    log_file_cf_update = os.path.join(logs_dir, "step2_classyfire_local_update.log")
    log_file_merge2 = os.path.join(logs_dir, "step2_merge_annotations_step2.log")
    log_file_filter = os.path.join(logs_dir, "step2_filter_after_merge.log") # log_file_filter 定义移到此处
    
    # 1. classyfire_add_batch.py
    py_classyfire_add_script = os.path.join(workflow_script_dir, "classyfire_add_batch.py")
    cmd_classyfire_add = ["python", py_classyfire_add_script]
    # log_file_cf_add 已在前面定义
    run_command(cmd_classyfire_add, mode_analysis_dir, f"Step 2 ({mode}) ClassyFire Add Batch: ", log_file_path=log_file_cf_add)

    # 2. classyfire_local_update.py
    py_classyfire_update_script = os.path.join(workflow_script_dir, "classyfire_local_update.py")
    cmd_classyfire_update = ["python", py_classyfire_update_script] 
    # log_file_cf_update 已在前面定义
    run_command(cmd_classyfire_update, mode_analysis_dir, f"Step 2 ({mode}) ClassyFire Local Update: ", log_file_path=log_file_cf_update)

    # 3. merge_multiple_annotations_step2.py
    py_merge_step2_script = os.path.join(workflow_script_dir, "merge_multiple_annotations_step2.py")
    cmd_merge_step2 = ["python", py_merge_step2_script] 
    # log_file_merge2 已在前面定义
    run_command(cmd_merge_step2, mode_analysis_dir, f"Step 2 ({mode}) Merge Annotations Step 2: ", log_file_path=log_file_merge2)

    # 4. filter_after_merge_annotation.py
    py_filter_script = os.path.join(workflow_script_dir, "filter_after_merge_annotation.py")
    input_filter_tsv = os.path.join(mode_analysis_dir, "merged_output_added_classyfire.tsv")
    output_clean_tsv = os.path.join(mode_analysis_dir, "merged_output_added_classyfire_clean.tsv")
    output_clean_filt_tsv = os.path.join(mode_analysis_dir, "merged_output_added_classyfire_clean_filt.tsv")
    cmd_filter = [
        "python", py_filter_script,
        input_filter_tsv, output_clean_tsv, output_clean_filt_tsv,
        str(num_samples)
    ]
    # log_file_filter 已在前面定义，此处不再重复定义
    run_command(cmd_filter, mode_analysis_dir, f"Step 2 ({mode}) Filter After Merge: ", log_file_path=log_file_filter)

    print(f"--- Finished Step 2 for {mode} mode ---")

def step3_combine(args):
    """
    执行步骤3：合并Pos和Neg结果，最终命名和过滤。
    """
    print("\n--- Starting Step 3 for Combine mode ---")
    project_analysis_dir = os.path.join(args.base_project_dir, args.project_name, "analysis")
    logs_dir = os.path.join(project_analysis_dir, "logs", "Combine") # Combine-specific log directory
    os.makedirs(logs_dir, exist_ok=True)
    print(f"Ensured log directory for Step 3 Combine: {logs_dir}")

    pos_analysis_dir = os.path.join(project_analysis_dir, "Pos")
    neg_analysis_dir = os.path.join(project_analysis_dir, "Neg")
    combine_analysis_dir = os.path.join(project_analysis_dir, "Combine")
    os.makedirs(combine_analysis_dir, exist_ok=True)

    workflow_script_dir = os.path.join(args.base_project_dir, "workflow_v2", "script")

    # 1. pos_neg_rbind.py
    py_rbind_script = os.path.join(workflow_script_dir, "pos_neg_rbind.py")
    pos_input_rbind = os.path.join(pos_analysis_dir, "merged_output_added_classyfire_clean_filt.tsv")
    neg_input_rbind = os.path.join(neg_analysis_dir, "merged_output_added_classyfire_clean_filt.tsv")
    combined_output_rbind = os.path.join(combine_analysis_dir, "combined_pos_neg.tsv")
    cmd_rbind = ["python", py_rbind_script, "-p", pos_input_rbind, "-n", neg_input_rbind, "-o", combined_output_rbind]
    log_file_rbind = os.path.join(logs_dir, "step3_pos_neg_rbind.log")
    run_command(cmd_rbind, combine_analysis_dir, "Step 3 (Combine) Rbind Pos/Neg: ", log_file_path=log_file_rbind)

    # 2. pass_filter_6s.py
    py_filter_6s_script = os.path.join(workflow_script_dir, "pass_filter_6s.py")
    output_6s = os.path.join(combine_analysis_dir, "combined_pos_neg_Passed.tsv")
    cmd_filter_6s = ["python", py_filter_6s_script, combined_output_rbind, output_6s]
    log_file_filter_6s = os.path.join(logs_dir, "step3_filter_6s.log")
    run_command(cmd_filter_6s, combine_analysis_dir, "Step 3 (Combine) Filter 6s bins: ", log_file_path=log_file_filter_6s)

    # 3. pass_filter_inchikey_prefix14.py
    py_filter_inchikey14_script = os.path.join(workflow_script_dir, "pass_filter_inchikey_prefix14.py")
    output_inchikey14 = os.path.join(combine_analysis_dir, "combined_pos_neg_Passed_inchikey14.tsv")
    cmd_filter_inchikey14 = ["python", py_filter_inchikey14_script, output_6s, output_inchikey14]
    log_file_filter_i14 = os.path.join(logs_dir, "step3_filter_inchikey14.log")
    run_command(cmd_filter_inchikey14, combine_analysis_dir, "Step 3 (Combine) Filter InChIKey Prefix 14: ", log_file_path=log_file_filter_i14)
    
    # 4. final_anno_name.py
    py_final_name_script = os.path.join(workflow_script_dir, "final_anno_name.py")
    cmd_final_name = ["python", py_final_name_script] 
    log_file_final_name = os.path.join(logs_dir, "step3_final_anno_name.log")
    run_command(cmd_final_name, combine_analysis_dir, "Step 3 (Combine) Final Annotation Name: ", log_file_path=log_file_final_name)

    # 5. simplify_table.py
    py_simplify_script = os.path.join(workflow_script_dir, "simplify_table.py")
    output_simplify_name = f"{args.project_name}-Metabolomics.txt"
    output_simplify = os.path.join(combine_analysis_dir, output_simplify_name)
    cmd_simplify = ["python", py_simplify_script, output_simplify]
    log_file_simplify = os.path.join(logs_dir, "step3_simplify_table.log")
    run_command(cmd_simplify, combine_analysis_dir, "Step 3 (Combine) Simplify Table: ", log_file_path=log_file_simplify)

    print(f"--- Finished Step 3 for Combine mode ---")
    print(f"\nWorkflow finished. Final simplified table: {output_simplify}")


def main():
    parser = argparse.ArgumentParser(description="Main workflow script for metabolomics data processing.")
    parser.add_argument("--project_name", required=True, help="Name of the project.")
    parser.add_argument("--base_project_dir", default="/mnt/NexData_file/home/zjh/project/Metabolomics_HILIC", help="Base directory for all projects.")
    parser.add_argument("--stage", type=int, choices=[1, 2], required=True, help="Workflow stage to execute (1: Step 1 only, 2: Step 2 and 3).")
    
    args = parser.parse_args()

    # Setup project base directories
    project_dir = os.path.join(args.base_project_dir, args.project_name)
    project_data_dir = os.path.join(project_dir, "data")
    project_analysis_dir = os.path.join(project_dir, "analysis")

    # Ensure base directories for data and analysis exist
    os.makedirs(project_data_dir, exist_ok=True)
    os.makedirs(project_analysis_dir, exist_ok=True)
    
    # Setup mode-specific subdirectories within data and analysis
    setup_directories(project_data_dir, ["Pos", "Neg"])
    setup_directories(project_analysis_dir, ["Pos", "Neg", "Combine", "logs", os.path.join("logs", "Pos"), os.path.join("logs", "Neg"), os.path.join("logs", "Combine")])

    print(f"Project: {args.project_name}")
    print(f"Base Directory: {args.base_project_dir}")
    print(f"Data Directory: {project_data_dir}")
    print(f"Analysis Directory: {project_analysis_dir}")
    print(f"Executing Stage: {args.stage}")

    if args.stage == 1:
        # Step 1: Run for Pos and Neg modes in parallel
        print("\n--- Starting Stage 1: Step 1 for Pos and Neg modes in parallel ---")
        # Informational check for mzXML files for stage 1
        pos_mzxml_dir_s1 = os.path.join(project_data_dir, "Pos")
        neg_mzxml_dir_s1 = os.path.join(project_data_dir, "Neg")
        num_samples_pos_s1 = len(glob.glob(os.path.join(pos_mzxml_dir_s1, "*.mzXML")))
        num_samples_neg_s1 = len(glob.glob(os.path.join(neg_mzxml_dir_s1, "*.mzXML")))
        print(f"Informational: Found {num_samples_pos_s1} .mzXML files in {pos_mzxml_dir_s1} for Pos mode (Stage 1).")
        print(f"Informational: Found {num_samples_neg_s1} .mzXML files in {neg_mzxml_dir_s1} for Neg mode (Stage 1).")
        if num_samples_pos_s1 == 0 or num_samples_neg_s1 == 0:
            print("Warning (Stage 1): No .mzXML files found in one or both mode directories. Proceeding with stage 1.")
        elif num_samples_pos_s1 != num_samples_neg_s1:
            print("Warning (Stage 1): Mismatch in the number of .mzXML files between Pos and Neg directories. Proceeding with stage 1.")

        process_pos = multiprocessing.Process(target=step1_mode, args=(args, "Pos"))
        process_neg = multiprocessing.Process(target=step1_mode, args=(args, "Neg"))

        process_pos.start()
        process_neg.start()

        process_pos.join()
        process_neg.join()
        print("\n--- Finished Stage 1: Step 1 for both Pos and Neg modes ---")
    
    elif args.stage == 2:
        print("\n--- Starting Stage 2: Step 2 and Step 3 ---")
        
        # Determine number of samples from .mzXML files - critical for stage 2
        pos_mzxml_dir = os.path.join(project_data_dir, "Pos")
        num_samples = len(glob.glob(os.path.join(pos_mzxml_dir, "*.mzXML")))

        print(f"Determined number of samples for Stage 2: {num_samples}")

        # Step 2: Run for Pos and Neg modes
        step2_mode(args, "Pos", num_samples)
        step2_mode(args, "Neg", num_samples)
        
        # Step 3: Combine results
        step3_combine(args)
        print("\n--- Finished Stage 2 ---")

    print("\nWorkflow processing complete for the selected stage.")

if __name__ == "__main__":
    main()
