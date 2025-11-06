use anyhow::Result;
use fastga_rs::{Config, FastGA};
/// Test that FastGA processes queries sequentially and completely
/// This is critical for our streaming API assumptions!

#[test]
fn test_fastga_processes_queries_sequentially() -> Result<()> {
    // Create test data with 3 distinct queries
    let temp_dir = tempfile::TempDir::new()?;

    // Create complex sequences that FastGA can align
    let queries_fasta = temp_dir.path().join("queries.fa");
    // Use yeast-like GC content and complexity
    let seq_a = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string()
        + "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        + "GATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAG"
        + "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
        + "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";
    let seq_b = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT".to_string()
        + "TATCGATGCATGCTAGCTACGATCGATGCTAGCTAGCATCGATCGATGCATGCTAGCTAG"
        + "CGATCGATGCTAGCTAGCATCGATCGATCGTAGCTAGCTGCATGCATGCTAGCTAGCTAG"
        + "ATGGCTAGCGATCGATGCATGCATCGATCGATGCTAGCTAGCATCGATGCTAGCTAGCAT"
        + "CGATCGATCGTAGCTAGCTAGCATGCATGCTAGCTAGCTAGCGATCGATCGATGCTAGCT";
    let seq_c = "GGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAA".to_string()
        + "TTCCAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAA"
        + "GCTTGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGG"
        + "ATCCCTCGAGCTCGAGCTCGAGCTCGAGCTCGAGCTCGAGCTCGAGCTCGAGCTCGAGCT"
        + "CGAGTCTAGATCTAGATCTAGATCTAGATCTAGATCTAGATCTAGATCTAGATCTAGATC";

    std::fs::write(
        &queries_fasta,
        format!(
            ">query_A\n{seq_a}\n>query_B\n{seq_b}\n>query_C\n{seq_c}\n"
        ),
    )?;

    // Target sequences with exact matches plus some variants
    let target_fasta = temp_dir.path().join("target.fa");

    // Create targets with good matches to each query
    let target_a = seq_a.clone(); // Exact match to query_A
    let target_b = seq_b.clone(); // Exact match to query_B
    let target_c = seq_c.clone(); // Exact match to query_C
                                  // Mixed target with partial matches from each
    let target_mixed = format!(
        "{}AAAAAA{}TTTTTT{}",
        &seq_a[..100],
        &seq_b[..100],
        &seq_c[..100]
    );

    std::fs::write(
        &target_fasta,
        format!(">target_1_matches_A\n{target_a}\n>target_2_matches_B\n{target_b}\n>target_3_matches_C\n{target_c}\n>target_4_mixed\n{target_mixed}\n")
    )?;

    // Use our fork-based FastGA API with single thread to ensure deterministic order
    let config = Config::builder()
        .num_threads(1)
        .min_alignment_length(20)
        .build();

    let aligner = FastGA::new(config)?;
    let alignments = aligner.align_files(&queries_fasta, &target_fasta)?;

    println!("FastGA output ({} alignments):", alignments.len());

    // Parse query names from alignments
    let mut query_order = Vec::new();
    let mut last_query = String::new();

    for alignment in &alignments.alignments {
        if alignment.query_name != last_query {
            query_order.push(alignment.query_name.clone());
            last_query = alignment.query_name.clone();
        }
    }

    // CRITICAL TEST: Queries must appear in order with no interleaving
    assert_eq!(
        query_order,
        vec!["query_A", "query_B", "query_C"],
        "FastGA must process queries sequentially!"
    );

    // Count alignments per query
    let mut query_a_count = 0;
    let mut query_b_count = 0;
    let mut query_c_count = 0;

    for alignment in &alignments.alignments {
        match alignment.query_name.as_str() {
            "query_A" => query_a_count += 1,
            "query_B" => query_b_count += 1,
            "query_C" => query_c_count += 1,
            _ => {}
        }
    }

    println!("\nAlignment counts:");
    println!("  query_A: {query_a_count} alignments");
    println!("  query_B: {query_b_count} alignments");
    println!("  query_C: {query_c_count} alignments");

    // Each query should have at least one alignment (to its matching target)
    assert!(query_a_count >= 1, "query_A should have alignments");
    assert!(query_b_count >= 1, "query_B should have alignments");
    assert!(query_c_count >= 1, "query_C should have alignments");

    // Verify no interleaving: Check that all query_A alignments come before query_B
    let mut seen_query_b = false;
    let mut seen_query_c = false;

    for alignment in &alignments.alignments {
        match alignment.query_name.as_str() {
            "query_A" => {
                assert!(
                    !seen_query_b && !seen_query_c,
                    "query_A alignment found after query_B or query_C!"
                );
            }
            "query_B" => {
                seen_query_b = true;
                assert!(!seen_query_c, "query_B alignment found after query_C!");
            }
            "query_C" => {
                seen_query_c = true;
            }
            _ => {}
        }
    }

    println!("\n✓ Query order verified: All queries processed sequentially");

    Ok(())
}

#[test]
fn test_query_completeness_with_multiple_targets() -> Result<()> {
    // Test that ALL alignments for a query are output before moving to next query
    let temp_dir = tempfile::TempDir::new()?;

    // Create complex sequences for better alignment
    let seq1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string()
        + "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        + "GATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAG"
        + "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
        + "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";
    let seq2 = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT".to_string()
        + "TATCGATGCATGCTAGCTACGATCGATGCTAGCTAGCATCGATCGATGCATGCTAGCTAG"
        + "CGATCGATGCTAGCTAGCATCGATCGATCGTAGCTAGCTGCATGCATGCTAGCTAGCTAG"
        + "ATGGCTAGCGATCGATGCATGCATCGATCGATGCTAGCTAGCATCGATGCTAGCTAGCAT"
        + "CGATCGATCGTAGCTAGCTAGCATGCATGCTAGCTAGCTAGCGATCGATCGATGCTAGCT";

    let queries_fasta = temp_dir.path().join("queries.fa");
    std::fs::write(
        &queries_fasta,
        format!(">query_1\n{seq1}\n>query_2\n{seq2}\n"),
    )?;

    // Multiple targets with good matches
    let target_fasta = temp_dir.path().join("target.fa");
    std::fs::write(
        &target_fasta,
        format!(
            ">target_A\n{}\n>target_B\n{}NNNN{}\n>target_C\n{}\n>target_D_no_match\n{}\n>target_E\n{}\n",
            seq1,  // Exact match to query_1
            &seq1[..150], &seq1[150..],  // Query_1 with gap
            &seq1[..200],  // Partial query_1
            ("GGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAA".to_string() +
             "TTCCAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAA" +
             "GCTTGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGG"),  // Different, no match
            seq2  // Exact match to query_2
        )
    )?;

    // Use our fork-based FastGA API
    let config = Config::builder()
        .num_threads(1)
        .min_alignment_length(10) // Lower threshold for test sequences
        .build();

    let aligner = FastGA::new(config)?;
    let alignments = aligner.align_files(&queries_fasta, &target_fasta)?;

    // Group alignments by query
    let mut current_query = String::new();
    let mut query_1_targets = Vec::new();
    let mut query_2_targets = Vec::new();

    for alignment in &alignments.alignments {
        let query = &alignment.query_name;
        let target = &alignment.target_name;

        if query == "query_1" {
            query_1_targets.push(target.clone());
            assert!(
                current_query.is_empty() || current_query == "query_1",
                "Query 1 alignments must not be interrupted"
            );
            current_query = "query_1".to_string();
        } else if query == "query_2" {
            query_2_targets.push(target.clone());
            assert!(
                current_query == "query_1" || current_query == "query_2",
                "Query 2 must come after Query 1"
            );
            current_query = "query_2".to_string();
        }
    }

    println!("\n✓ Query completeness verified:");
    println!("  Query 1 aligned to: {query_1_targets:?}");
    println!("  Query 2 aligned to: {query_2_targets:?}");

    if !query_1_targets.is_empty() && !query_2_targets.is_empty() {
        println!("  All query_1 alignments completed before query_2 started");
    }

    Ok(())
}

#[test]
fn test_parallel_threads_still_sequential_queries() -> Result<()> {
    // Even with multiple threads, queries should be processed sequentially
    let temp_dir = tempfile::TempDir::new()?;

    // Create simple test data
    let queries_fasta = temp_dir.path().join("queries.fa");
    std::fs::write(
        &queries_fasta,
        ">query_X\nACGTACGTACGTACGT\n>query_Y\nTGCATGCATGCATGCA\n",
    )?;

    let target_fasta = temp_dir.path().join("target.fa");
    std::fs::write(
        &target_fasta,
        ">target_1\nACGTACGTACGTACGT\n>target_2\nTGCATGCATGCATGCA\n",
    )?;

    // Test with multiple threads
    let config = Config::builder()
        .num_threads(4) // Multiple threads
        .min_alignment_length(10)
        .build();

    let aligner = FastGA::new(config)?;
    let alignments = aligner.align_files(&queries_fasta, &target_fasta)?;

    // Verify query order even with parallel processing
    let mut last_query = String::new();
    for alignment in &alignments.alignments {
        if alignment.query_name != last_query {
            println!("Processing query: {}", alignment.query_name);
            assert!(
                last_query.is_empty() || alignment.query_name > last_query,
                "Queries must be processed in order even with multiple threads"
            );
            last_query = alignment.query_name.clone();
        }
    }

    Ok(())
}
