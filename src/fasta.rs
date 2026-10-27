//! Record-preserving DNA FASTA input for exact matching.

use std::collections::HashSet;
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: Vec<u8>,
}

/// Parse FASTA without joining records or changing their biological coordinates.
/// ASCII whitespace is ignored. DNA/RNA IUPAC letters are retained in uppercase;
/// the matching layer treats every non-ACGT letter as a nonmatching barrier.
/// Empty records, duplicate IDs, missing headers and other characters are errors.
pub fn parse_fasta_records(text: &str) -> Result<Vec<FastaRecord>, String> {
    let mut records: Vec<FastaRecord> = Vec::new();
    let mut ids = HashSet::new();
    for (line_index, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if records.last().is_some_and(|r| r.sequence.is_empty()) {
                return Err(format!("Empty FASTA record before line {}", line_index + 1));
            }
            let id = header
                .split_whitespace()
                .next()
                .ok_or_else(|| format!("Missing FASTA identifier on line {}", line_index + 1))?;
            if !ids.insert(id.to_owned()) {
                return Err(format!("Duplicate FASTA identifier: {id}"));
            }
            records.push(FastaRecord {
                id: id.to_owned(),
                sequence: Vec::new(),
            });
        } else {
            let record = records.last_mut().ok_or_else(|| {
                format!(
                    "Sequence before first FASTA header on line {}",
                    line_index + 1
                )
            })?;
            for base in line.bytes().filter(|b| !b.is_ascii_whitespace()) {
                let base = base.to_ascii_uppercase();
                if !b"ACGTRYSWKMBDHVNU".contains(&base) {
                    return Err(format!(
                        "Invalid DNA FASTA byte on line {}: {base}",
                        line_index + 1
                    ));
                }
                record.sequence.push(base);
            }
        }
    }
    if records.is_empty() {
        return Err("FASTA file contains no records".to_owned());
    }
    if records.last().is_some_and(|r| r.sequence.is_empty()) {
        return Err("Last FASTA record is empty".to_owned());
    }
    Ok(records)
}

pub fn read_fasta_records(path: impl AsRef<Path>) -> Result<Vec<FastaRecord>, String> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    parse_fasta_records(&text).map_err(|e| format!("{}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn records_case_wrapping_whitespace_and_crlf() {
        let records =
            parse_fasta_records(">r1 description\r\nac gt\r\nNN\r\n\r\n>r2\r\ntg\tca\r\n").unwrap();
        assert_eq!(
            records,
            vec![
                FastaRecord {
                    id: "r1".into(),
                    sequence: b"ACGTNN".to_vec()
                },
                FastaRecord {
                    id: "r2".into(),
                    sequence: b"TGCA".to_vec()
                },
            ]
        );
    }

    #[test]
    fn invalid_input_is_not_silently_converted_to_n() {
        for input in [
            "",
            "ACGT",
            ">\nACGT",
            ">a\n",
            ">a\n>b\nA",
            ">a\nA\n>a\nT",
            ">a\nAC?T",
            ">a\nAC\0T",
        ] {
            assert!(parse_fasta_records(input).is_err(), "accepted {input:?}");
        }
    }
}
