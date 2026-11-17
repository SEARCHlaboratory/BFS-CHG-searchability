from CHGGenerator import chg_generator  
import time
import re
import random 
import numpy as np
from itertools import product
import statistics
import csv

# random.seed(9 + run)
# np.random.seed(1 + run)
#vars

num_runs = 30
sizes = [50, 100, 200]
connectedness_levels = [0.2, 0.6, 1.0]

factorial_design = {
    "Simple Full Acyclic": {
        "multiedges": [0],
        "partiality": [0],
        "cyclicity": [0]
    },
    "Simple Partial Acyclic": {
        "multiedges": [0],
        "partiality": [0.5, 0.8],
        "cyclicity": [0]
    },
    "Simple Partial Cyclic": {
        "multiedges": [0],
        "partiality": [0.5, 0.8],
        "cyclicity": [0.5, 0.8]
    },
    "Non-simple Full Acyclic": {
        "multiedges": [0.5, 0.8],
        "partiality": [0],
        "cyclicity": [0]
    },
    "Non-simple Partial Acyclic": {
        "multiedges": [0.5, 0.8],
        "partiality": [0.5, 0.8],
        "cyclicity": [0]
    },
    "Non-simple Partial Cyclic": {
        "multiedges": [0.5, 0.8],
        "partiality": [0.5, 0.8],
        "cyclicity": [0.5, 0.8]
    }
}

cases = []
all_results = []
case_averages = []


for structure, factors in factorial_design.items():
    for multi, partial, cyclic, size, connected in product(
        factors["multiedges"],
        factors["partiality"],
        factors["cyclicity"],
        sizes,
        connectedness_levels
    ):
        cases.append({
            "structure": structure,
            "num_nodes": size,
            "connectedness": connected,
            "prop_multiedges": multi,
            "prop_partial_edges": partial,
            "prop_cycles": cyclic
        })

print("Total Combination of Cases:", len(cases))
print("Total Planned Runs:", len(cases) * num_runs)


inputs = {"AA": 5, "AB": 7}

solved_chgs = 0
formed_chgs = 0
tests = 0



for case_number, case in enumerate(cases, start=1):

    num_nodes = case["num_nodes"]
    connectedness = case["connectedness"]
    prop_partial_edges = case["prop_partial_edges"]
    prop_cycles = case["prop_cycles"]
    prop_multiedges = case["prop_multiedges"]

#========================================= print
    
#===========================================================
    solved_chgs = 0
    formed_chgs = 0
    tests = 0

    results = []
    logged_times = []
    unlogged_times = []

    print("\n" + "=" * 65)
    print(f"COMB. CASE {case_number}/{len(cases)}")
    print("CHG Structure:", case["structure"])
    print("Nodes:", num_nodes)
    print("Connectedness:", connectedness)
    print("Partiality:", prop_partial_edges)
    print("Cyclicity:", prop_cycles)
    print("Multiedges:", prop_multiedges)
    print("=" * 65)

    for run in range(num_runs):
        random.seed(9 + run)
        np.random.seed(1 + run)
    # print(f"\nStarting CHG {run+1}/{num_runs}", flush=True)
    # Generate a new CHG for each run
        try:
            # print("Generating CHG...", flush=True)
            hg = chg_generator(num_nodes,connectedness,prop_partial_edges,prop_cycles,prop_multiedges)
            formed_chgs += 1
        except Exception as e:
            print(f"\nRUN {run+1} - Generation Failed:", e)
            continue
    # print("CHG generated. Starting logged BFS...", flush=True)

    # Select target from generated CHG
        target = random.choice(list(hg.nodes.keys()))

        total_nodes = len(hg.nodes)
        total_edges = len(hg.edges)
        connected_nodes = len(hg.targeted_nodes)
        num_edges = int(num_nodes * connectedness)
        density = (len(hg.edges) / total_nodes) 

        # Clear log
        with open("constrainthg.log", "w") as file:
            file.write("")

    #logged run
        start = time.perf_counter()
        try:
            t2 = hg.solve(target, inputs, search_depth=num_edges*3, logging_level=9)
            logged_times.append(time.perf_counter() - start)
        except Exception as e:
            t2 = None
            status = str(e)
        runtime = time.perf_counter() - start
        
        # print("Logged BFS finished. Starting unlogged BFS...", flush=True)

        with open("constrainthg.log", "r") as file:
                log = file.read()

        # Read log before running without logging
        NLstart = time.perf_counter()
        try:
            t = hg.solve(target, inputs, search_depth=num_edges*3,)
            unlogged_times.append(time.perf_counter() - NLstart)
            status = "Solved" if t is not None else "Not reached"
        except Exception as e:
            t = None
            status = str(e)
        unlogged_runtime = time.perf_counter() - NLstart

        tests += 1

        if t is not None:
            solved_chgs += 1
            status = "Solved"
        else:
            status = "Not reached"


        section = log.split("Explored edges(")[-1]
        edge_matches = re.findall(r"<(.+?)>:\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)",section)

        explored_edges = sum(1 for e in edge_matches if int(e[1]) > 0)
        processed_edges = sum(1 for e in edge_matches if int(e[2]) > 0)
        valid_edges = sum(1 for e in edge_matches if int(e[3]) > 0)

        total_explored = sum(int(x[1]) for x in edge_matches)
        total_processed = sum(int(x[2]) for x in edge_matches)
        total_valid_solutions = sum(int(x[3]) for x in edge_matches)

        # Extract visited nodes and states
        visited = re.findall(r"Exploring <([^>]+)>", log)
        unique_nodes = list(dict.fromkeys(state.split("#")[0] for state in visited))
        total_states = len(visited)

        # Calculate percentages
        percent_explored = ((explored_edges / total_edges) * 100 if total_edges > 0 else 0)
        percent_processed = ((processed_edges / explored_edges) * 100 if explored_edges > 0 else 0)

        percent_valid = ((valid_edges / processed_edges) * 100 if processed_edges > 0 else 0)
        s_path_len = len(t.values)-1 if t is not None else None

        run_result = {
            "Case": case_number,
            "Run": run + 1,
            "Structure": case["structure"],
            "Node Size": num_nodes,
            "Connectedness": connectedness,
            "Partiality": prop_partial_edges,
            "Cyclicity": prop_cycles,
            "Multiedges": prop_multiedges,
            "Search Status": status,
            "Total Nodes": total_nodes,
            "Connected Nodes": connected_nodes,
            "Total Edges": total_edges,
            "Density": density,
            "Unique Nodes": len(unique_nodes),
            "States Visited": total_states,
            "Total Explored": total_explored if edge_matches else None,
            "Total Processed": total_processed if edge_matches else None,
            "Total Valid Solutions": total_valid_solutions if edge_matches else None,
            "Percent Explored": percent_explored,
            "Percent Processed": percent_processed,
            "Percent Valid": percent_valid,
            "Target Node": target,
            "Target Value": t.value if t is not None else None,
            "Solution Path Length": s_path_len,
            "Logged Runtime": runtime,
            "Unlogged Runtime": unlogged_runtime
        }

    #store individual run results
        
        results.append({
            "runtime": runtime,
            "unlogged_runtime": unlogged_runtime,
            "status": status,
            "total_nodes": total_nodes,
            "connected_nodes": len(hg.targeted_nodes),
            "density": density,
            "total_edges": total_edges,
            "unique_nodes": len(unique_nodes),
            "total_states": total_states,
            "target_node": target,
            "percent_explored": percent_explored,
            "percent_processed": percent_processed,
            "percent_valid": percent_valid,
            "total_explored": total_explored if edge_matches else None,
            "total_processed": total_processed if edge_matches else None,
            "total_valid_solutions": total_valid_solutions if edge_matches else None,
            "solution_path_length": s_path_len
        })
        all_results.append(run_result) #for printer
#===================================================
        
#===================================================
        

    #print individual run results

        print(f"\nRUN {run+1}")
        print("Search Status:", status)
        print("Total Nodes:", total_nodes)
        print("Connected Nodes:", connected_nodes)
        print("Total Edges:", total_edges)
        print("Density:", density)
        print("Unique Nodes:", unique_nodes)
        print("Total Unique Nodes:", len(unique_nodes))
        print("Total States Visited:", total_states)
        print("Total Explored:", total_explored if edge_matches else "N/A")
        print("Total Processed:", total_processed if edge_matches else "N/A")
        print("Total Valid Solutions:", total_valid_solutions if edge_matches else "N/A")
        print(f"Percent Explored: {percent_explored:.2f}%")
        print(f"Percent Processed: {percent_processed:.2f}%")
        print(f"Percent Valid: {percent_valid:.2f}%")
        print("Target Node:", target)
        print("Target Value:", t.value if t is not None else "Not reached")
        print("Solution Path Length:", s_path_len if s_path_len is not None else "Not reached")
        print("LRuntime:", runtime)
        print("NLRuntime:", unlogged_runtime)

        if t is not None:
           print("\nSOLUTION PATH VALUES")
           for node, values in reversed(list(t.values.items())):
               print(f"{node}: {values}")
#calculate averages

    def average(metric):
        values = [r[metric] for r in results if r[metric] is not None]
        return statistics.mean(values) if values else None

    avg_Lruntime = average("runtime")
    avg_NLruntime = average("unlogged_runtime")
    avg_density = average("density")
    avg_connected = average("connected_nodes")
    avg_unique = average("unique_nodes")
    avg_states = average("total_states")
    avg_percent = average("percent_explored")
    avg_per_processed = average("percent_processed")
    avg_per_valid = average("percent_valid")
    avg_explored = average("total_explored")
    avg_processed = average("total_processed")
    avg_valid = average("total_valid_solutions")
    avg_path_len = average("solution_path_length")
    avg_nodes = average("total_nodes")
    avg_edges = average("total_edges")
    # avg_req_edges = average("num_edges")

   #Printer
    case_averages.append({
        "Case": case_number,
        "Structure": case["structure"],
        "Node Size": num_nodes,
        "Connectedness": connectedness,
        "Partiality": prop_partial_edges,
        "Cyclicity": prop_cycles,
        "Multiedges": prop_multiedges,
        "Runs Requested": num_runs,
        "Created CHGs": formed_chgs,
        "Successful Runs": solved_chgs,
        "Average Total Nodes": avg_nodes,
        "Average Connected Nodes": avg_connected,
        "Average Total Edges": avg_edges,
        "Average Density": avg_density,
        "Average Unique Nodes": avg_unique,
        "Average States Visited": avg_states,
        "Average Total Explored": avg_explored,
        "Average Total Processed": avg_processed,
        "Average Total Valid Solutions": avg_valid,
        "Average Percent Explored": avg_percent,
        "Average Percent Processed": avg_per_processed,
        "Average Percent Valid": avg_per_valid,
        "Average Path Length": avg_path_len,
        "Average Logged Runtime": avg_Lruntime,
        "Average Unlogged Runtime": avg_NLruntime,
        "Average Logging Overhead": (
            avg_Lruntime - avg_NLruntime
            if avg_Lruntime is not None and avg_NLruntime is not None
            else None
        )
    })

#==========================================
   
#===========================================


#print overall avg results

    print("\n" + "=" * 65)
    print(f"AVERAGE RESULTS ACROSS {len(results)} CHGS")
    print("=" * 65)

    print("Total Runs Requested:", num_runs)
    print("Created CHGs:", formed_chgs)
    print("Successful Runs:", sum(r["status"] == "Solved" for r in results))
    print("Proportion Partial Edges:", prop_partial_edges)
    print("Proportion Cycles:", prop_cycles)
    print("Proportion Multiedges:", prop_multiedges)
    print("Average Total Nodes:", f"{avg_nodes:.2f}" if avg_nodes is not None else "N/A")
    print("Average Connected Nodes:", f"{avg_connected:.2f}" if avg_connected is not None else "N/A")
    print("Average Total Edges:", f"{avg_edges:.2f}" if avg_edges is not None else "N/A")
    print("Average Density:", f"{avg_density:.2f}" if avg_density is not None else "N/A")
    print("Average Unique Nodes:", f"{avg_unique:.2f}" if avg_unique is not None else "N/A")
    print("Average States Visited:", f"{avg_states:.2f}" if avg_states is not None else "N/A")
    print("Average Total Explored:", f"{avg_explored:.2f}" if avg_explored is not None else "N/A")
    print("Average Total Processed:", f"{avg_processed:.2f}" if avg_processed is not None else "N/A")
    print("Average Total Valid Solutions:", f"{avg_valid:.2f}" if avg_valid is not None else "N/A")
    print("Average Percent Explored:", f"{avg_percent:.2f}%" if avg_percent is not None else "N/A")
    print("Average Percent Processed:", f"{avg_per_processed:.2f}%" if avg_per_processed is not None else "N/A")
    print("Average Percent Valid:", f"{avg_per_valid:.2f}%" if avg_per_valid is not None else "N/A")
    print("Average Path Length:", f"{avg_path_len:.2f}" if avg_path_len is not None else "N/A")
    print("Average Logged Runtime:", avg_Lruntime)
    print("Average Unlogged Runtime:", avg_NLruntime)
    print(
        "Average Logging Overhead:",
        avg_Lruntime - avg_NLruntime
        if avg_Lruntime is not None and avg_NLruntime is not None
        else "N/A"
    )
    print(f"Created {formed_chgs} CHGS.\nSolved {solved_chgs} CHGs out of {tests} tests")

#print individual results

if all_results:
    with open("CHG_Individual_Results.csv", "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=all_results[0].keys()
        )
        writer.writeheader()
        writer.writerows(all_results)

    print("Saved: CHG_Individual_Results_Broad.csv")


#export avg results

if case_averages:
    with open("CHG_Average_Results_Broad.csv", "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=case_averages[0].keys()
        )
        writer.writeheader()
        writer.writerows(case_averages)

    print("Saved: CHG_Average_Results.csv")
