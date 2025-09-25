use anyhow::Result;
use fastga_rs::{Config, FastGA};
/// Test that FastGA processes queries sequentially and completely
/// This is critical for our streaming API assumptions!
use std::process::Command;

#[test]
fn test_fastga_processes_queries_sequentially() -> Result<()> {
    // Create test data with 3 distinct queries
    let temp_dir = tempfile::TempDir::new()?;

    // Queries with unique sequences that will align differently
    let queries_fasta = temp_dir.path().join("queries.fa");
    std::fs::write(
        &queries_fasta,
        ">query_A\n\
         ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
         >query_B\n\
         TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n\
         >query_C\n\
         AAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGG\n",
    )?;

    // Target with sequences that match each query differently
    let target_fasta = temp_dir.path().join("target.fa");
    std::fs::write(
        &target_fasta,
        ">target_1_matches_A\n\
         ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
         >target_2_matches_B\n\
         TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n\
         >target_3_matches_C\n\
         AAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGG\n\
         >target_4_mixed\n\
         ACGTACGTTGCATGCAAAAATTTTCCCCGGGGACGTACGTTGCATGCAAAAATTTTCCCCGGGGACGTACGTTGCATGCA\n",
    )?;

    // Convert to GDB - FAtoGDB syntax is: FAtoGDB <input.fa> <output_name>
    // It will create output_name.gdb
    let query_gdb_base = temp_dir.path().join("queries");
    let target_gdb_base = temp_dir.path().join("target");

    let fa_output = Command::new("deps/fastga/FAtoGDB")
        .arg(&queries_fasta)
        .arg(&query_gdb_base)
        .output()?;

    if !fa_output.status.success() {
        eprintln!(
            "FAtoGDB failed for queries: {}",
            String::from_utf8_lossy(&fa_output.stderr)
        );
    }

    let fa_output2 = Command::new("deps/fastga/FAtoGDB")
        .arg(&target_fasta)
        .arg(&target_gdb_base)
        .output()?;

    if !fa_output2.status.success() {
        eprintln!(
            "FAtoGDB failed for target: {}",
            String::from_utf8_lossy(&fa_output2.stderr)
        );
    }

    // List what files were created
    println!("Files in temp dir:");
    for entry in std::fs::read_dir(temp_dir.path())? {
        let entry = entry?;
        println!("  {:?}", entry.path());
    }

    // Run FastGA
    let output = Command::new("deps/fastga/FastGA")
        .arg("-T1") // Single thread to ensure deterministic order
        .arg("-pafx") // PAF with extended CIGAR
        .arg(&query_gdb_base) // Will look for queries.gdb
        .arg(&target_gdb_base) // Will look for target.gdb
        .output()?;

    // Debug output
    if !output.status.success() {
        eprintln!("FastGA failed: {}", String::from_utf8_lossy(&output.stderr));
    }

    let paf = String::from_utf8(output.stdout)?;
    println!("FastGA output ({} bytes):", paf.len());
    println!("{}", paf);
    let lines: Vec<&str> = paf.lines().collect();

    // Parse query names from PAF output
    let mut query_order = Vec::new();
    let mut last_query = String::new();

    for line in &lines {
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let query_name = fields[0];
        if query_name != last_query {
            query_order.push(query_name.to_string());
            last_query = query_name.to_string();
        }
    }

    // CRITICAL TEST: Queries must appear in order with no interleaving
    assert_eq!(
        query_order,
        vec!["query_A", "query_B", "query_C"],
        "FastGA must process queries sequentially!"
    );

    // Count alignments per query
    let mut query_a_lines = Vec::new();
    let mut query_b_lines = Vec::new();
    let mut query_c_lines = Vec::new();

    for (i, line) in lines.iter().enumerate() {
        if line.starts_with("query_A") {
            query_a_lines.push(i);
        } else if line.starts_with("query_B") {
            query_b_lines.push(i);
        } else if line.starts_with("query_C") {
            query_c_lines.push(i);
        }
    }

    // Verify no interleaving: All query_A lines come before query_B, etc.
    if !query_a_lines.is_empty() && !query_b_lines.is_empty() {
        assert!(
            *query_a_lines.last().unwrap() < *query_b_lines.first().unwrap(),
            "Query A alignments must complete before Query B starts"
        );
    }
    if !query_b_lines.is_empty() && !query_c_lines.is_empty() {
        assert!(
            *query_b_lines.last().unwrap() < *query_c_lines.first().unwrap(),
            "Query B alignments must complete before Query C starts"
        );
    }

    println!("\n✓ FastGA processes queries sequentially:");
    println!("  - Query A: lines {:?}", query_a_lines);
    println!("  - Query B: lines {:?}", query_b_lines);
    println!("  - Query C: lines {:?}", query_c_lines);

    Ok(())
}

#[test]
fn test_query_completeness_with_multiple_targets() -> Result<()> {
    // Test that ALL alignments for a query are output before moving to next query
    let temp_dir = tempfile::TempDir::new()?;

    // One query that should match multiple targets
    let queries_fasta = temp_dir.path().join("queries.fa");
    std::fs::write(
        &queries_fasta,
        ">query_1\n\
         ACGTACGTACGTACGT\n\
         >query_2\n\
         TGCATGCATGCATGCA\n",
    )?;

    // Multiple targets with partial matches to query_1
    let target_fasta = temp_dir.path().join("target.fa");
    std::fs::write(
        &target_fasta,
        ">target_A\n\
         ACGTACGTACGTACGT\n\
         >target_B\n\
         NNNNACGTACGTACGTACGTNNNN\n\
         >target_C\n\
         ACGTACGTNNNNNNNNACGTACGT\n\
         >target_D_no_match\n\
         TTTTTTTTTTTTTTTTTTTTTTTT\n\
         >target_E\n\
         TGCATGCATGCATGCA\n",
    )?;

    // Convert and run
    let query_gdb_base = temp_dir.path().join("queries");
    let target_gdb_base = temp_dir.path().join("target");

    Command::new("deps/fastga/FAtoGDB")
        .arg(&queries_fasta)
        .arg(&query_gdb_base)
        .output()?;

    Command::new("deps/fastga/FAtoGDB")
        .arg(&target_fasta)
        .arg(&target_gdb_base)
        .output()?;

    let output = Command::new("deps/fastga/FastGA")
        .arg("-T1")
        .arg("-pafx")
        .arg("-l50") // Lower min length for test
        .arg(&query_gdb_base)
        .arg(&target_gdb_base)
        .output()?;

    let paf = String::from_utf8(output.stdout)?;

    // Group alignments by query
    let mut current_query = String::new();
    let mut query_1_targets = Vec::new();
    let mut query_2_targets = Vec::new();

    for line in paf.lines() {
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let query = fields[0];
        let target = fields[5];

        if query == "query_1" {
            query_1_targets.push(target);
            assert!(
                current_query.is_empty() || current_query == "query_1",
                "Query 1 alignments must not be interrupted"
            );
            current_query = "query_1".to_string();
        } else if query == "query_2" {
            query_2_targets.push(target);
            assert!(
                current_query == "query_1" || current_query == "query_2",
                "Query 2 must come after Query 1"
            );
            current_query = "query_2".to_string();
        }
    }

    println!("\n✓ Query completeness verified:");
    println!("  Query 1 aligned to: {:?}", query_1_targets);
    println!("  Query 2 aligned to: {:?}", query_2_targets);
    println!("  All query_1 alignments completed before query_2 started");

    Ok(())
}

#[test]
#[ignore = "Requires FastGA binaries to be built"]
fn test_parallel_threads_still_sequential_queries() -> Result<()> {
    // Even with multiple threads, queries should be processed sequentially
    // (threads parallelize within a query, not across queries)

    let temp_dir = tempfile::TempDir::new()?;

    // Create queries
    let queries_fasta = temp_dir.path().join("queries.fa");
    let mut query_content = String::new();
    for i in 0..10 {
        query_content.push_str(&format!(">query_{:02}\n", i));
        query_content.push_str(&"ACGT".repeat(25));
        query_content.push('\n');
    }
    std::fs::write(&queries_fasta, query_content)?;

    // Run with 8 threads
    let output = Command::new("deps/fastga/FastGA")
        .arg("-T8") // Multiple threads
        .arg("-pafx")
        .arg(&queries_fasta)
        .arg(&queries_fasta) // Self-alignment
        .output()?;

    let paf = String::from_utf8(output.stdout)?;

    // Verify query order is maintained
    let mut last_query_num = -1;
    for line in paf.lines() {
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 1 {
            continue;
        }

        if let Some(num_str) = fields[0].strip_prefix("query_") {
            let num: i32 = num_str.parse()?;
            assert!(
                num >= last_query_num,
                "Query {} appeared after query {}",
                num,
                last_query_num
            );
            last_query_num = num;
        }
    }

    println!("\n✓ Multi-threaded FastGA still processes queries in order");

    Ok(())
}
