use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct GenomicStats {
    pub num_sequences: usize,
    pub total_length: usize,
    pub mean_length: f64,
    pub n50: usize,
    pub n90: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub gc_content: f64,
}

impl GenomicStats {
    pub fn new(sequences: &[Vec<u8>]) -> Self {
        if sequences.is_empty() {
            return GenomicStats {
                num_sequences: 0,
                total_length: 0,
                mean_length: 0.0,
                n50: 0,
                n90: 0,
                min_length: 0,
                max_length: 0,
                gc_content: 0.0,
            };
        }

        let lengths: Vec<usize> = sequences.iter().map(|seq| seq.len()).collect();
        let total_length: usize = lengths.iter().sum();
        let num_sequences = sequences.len();
        let mean_length = total_length as f64 / num_sequences as f64;
        
        let min_length = *lengths.iter().min().unwrap_or(&0);
        let max_length = *lengths.iter().max().unwrap_or(&0);

        // Calculate N50 and N90
        let mut sorted_lengths = lengths.clone();
        sorted_lengths.sort_by(|a, b| b.cmp(a)); // Sort in descending order
        
        let n50 = Self::calculate_nx(&sorted_lengths, 50.0);
        let n90 = Self::calculate_nx(&sorted_lengths, 90.0);

        // Calculate GC content
        let gc_count = sequences.iter()
            .flat_map(|seq| seq.iter())
            .filter(|&&base| base == b'G' || base == b'g' || base == b'C' || base == b'c')
            .count();
        let gc_content = if total_length > 0 {
            gc_count as f64 / total_length as f64 * 100.0
        } else {
            0.0
        };

        GenomicStats {
            num_sequences,
            total_length,
            mean_length,
            n50,
            n90,
            min_length,
            max_length,
            gc_content,
        }
    }

    fn calculate_nx(lengths: &[usize], percentage: f64) -> usize {
        let total_length: usize = lengths.iter().sum();
        let target_length = (total_length as f64 * percentage / 100.0).round() as usize;
        
        let mut current_length = 0;
        for &len in lengths {
            current_length += len;
            if current_length >= target_length {
                return len;
            }
        }
        
        lengths.first().copied().unwrap_or(0)
    }

    pub fn print_stats(&self, label: &str) {
        println!("{} Statistics:", label);
        println!("  Number of sequences: {}", self.num_sequences);
        println!("  Total length: {}", self.total_length);
        println!("  Mean length: {:.2}", self.mean_length);
        println!("  Min length: {}", self.min_length);
        println!("  Max length: {}", self.max_length);
        println!("  N50: {}", self.n50);
        println!("  N90: {}", self.n90);
        println!("  GC content: {:.2}%", self.gc_content);
        println!();
    }
}

pub fn parse_fasta(filename: &str) -> Vec<Vec<u8>> {
    let content = std::fs::read_to_string(filename)
        .expect("Could not read file");
    
    let mut sequences = Vec::new();
    let mut current_seq = Vec::new();
    
    for line in content.lines() {
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                sequences.push(current_seq);
                current_seq = Vec::new();
            }
        } else {
            current_seq.extend_from_slice(line.as_bytes());
        }
    }
    
    // Don't forget the last sequence
    if !current_seq.is_empty() {
        sequences.push(current_seq);
    }
    
    // Convert to uppercase and validate DNA sequence
    for seq in &mut sequences {
        for base in seq {
            *base = match *base {
                b'a' | b'A' => b'A',
                b'c' | b'C' => b'C',
                b'g' | b'G' => b'G',
                b't' | b'T' => b'T',
                b'n' | b'N' => b'N',
                _ => b'N', // Default to N for non-standard bases
            };
        }
    }
    
    sequences
}
