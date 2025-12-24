use crate::Match;

#[derive(Debug, Clone)]
pub enum OutputFormat {
    Default,
    Delta,
    Paf,
    Sam,
}

impl OutputFormat {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "delta" => Some(OutputFormat::Delta),
            "paf" => Some(OutputFormat::Paf),
            "sam" => Some(OutputFormat::Sam),
            _ => None,
        }
    }
}

pub fn print_matches_in_format(matches: &[Match], query_file: &str, format: &OutputFormat, reference_seq: &[u8], query_seq: &[u8]) {
    match format {
        OutputFormat::Default => print_matches_default(matches, query_file),
        OutputFormat::Delta => print_matches_delta(matches, query_file, reference_seq, query_seq),
        OutputFormat::Paf => print_matches_paf(matches, query_file, reference_seq, query_seq),
        OutputFormat::Sam => print_matches_sam(matches, query_file, reference_seq, query_seq),
    }
}

fn print_matches_default(matches: &[Match], query_file: &str) {
    println!("> Query: {}", query_file);
    for m in matches {
        println!("  Ref: {}  Query: {}  Len: {}", m.ref_pos + 1, m.query_pos + 1, m.len);
    }
}

fn print_matches_delta(matches: &[Match], query_file: &str, reference_seq: &[u8], _query_seq: &[u8]) {
    // Print header for delta format
    println!("NUCMER");
    println!("NUCMER");
    
    // Print reference and query file names
    println!("> Reference Query");
    
    for m in matches {
        // Delta format: ref_start ref_end query_start query_end ref_len query_len match_len
        let ref_start = m.ref_pos + 1;  // 1-based indexing
        let ref_end = m.ref_pos + m.len;
        let query_start = m.query_pos + 1;  // 1-based indexing
        let query_end = m.query_pos + m.len;
        
        let ref_len = reference_seq.len();
        let query_len = _query_seq.len();
        
        println!("{} {} {} {} {} {} {}", 
                 ref_start, ref_end, query_start, query_end, ref_len, query_len, m.len);
    }
}

fn print_matches_paf(matches: &[Match], query_file: &str, reference_seq: &[u8], query_seq: &[u8]) {
    for m in matches {
        // PAF format: query_name, query_length, query_start, query_end, 
        // strand, ref_name, ref_length, ref_start, ref_end, 
        // matching_bases, alignment_length, mapping_quality
        
        let query_name = query_file;
        let query_length = query_seq.len();
        let query_start = m.query_pos;
        let query_end = m.query_pos + m.len;
        
        let strand = "+"; // For simplicity, assuming forward strand
        
        let ref_name = "reference"; // Using a generic name
        let ref_length = reference_seq.len();
        let ref_start = m.ref_pos;
        let ref_end = m.ref_pos + m.len;
        
        let matching_bases = m.len; // Assuming all bases match for simplicity
        let alignment_length = m.len;
        let mapping_quality = 60; // Default mapping quality
        
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                 query_name, query_length, query_start, query_end,
                 strand, ref_name, ref_length, ref_start, ref_end,
                 matching_bases, alignment_length, mapping_quality);
    }
}

fn print_matches_sam(matches: &[Match], query_file: &str, reference_seq: &[u8], query_seq: &[u8]) {
    // Print SAM header if this is the first output
    println!("@HD\tVN:1.6");
    println!("@SQ\tSN:reference\tLN:{}", reference_seq.len());
    
    for m in matches {
        // SAM format: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
        
        let qname = query_file; // Query template NAME
        let flag = 0; // bitwise FLAG (0 for forward strand, unmated)
        let rname = "reference"; // Reference sequence NAME
        let pos = m.ref_pos + 1; // 1-based leftmost mapping POSition
        let mapq = 60; // MAPping Quality
        let cigar = format!("{}M", m.len); // CIGAR string
        let rnext = "*"; // Ref. name of the mate/next read
        let pnext = 0; // Position of the mate/next read
        let tlen = 0; // observed Template LENgth
        let seq = String::from_utf8_lossy(&query_seq[m.query_pos..m.query_pos + m.len]); // segment SEQuence
        let qual = "*"; // ASCII of Phred-scaled base QUALity+33
        
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                 qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual);
    }
}
