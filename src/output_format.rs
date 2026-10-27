//! Output for contig-aware, strand-aware exact matches.
use crate::{FastaRecord, RecordMatch};
use std::io::{self, Write};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Default,
    Matches,
    Paf,
    Delta,
    Sam,
}

impl OutputFormat {
    pub fn from_str(value: &str) -> Option<Self> {
        match value.to_ascii_lowercase().as_str() {
            "default" => Some(Self::Default),
            "matches" | "tsv" => Some(Self::Matches),
            "paf" => Some(Self::Paf),
            "delta" => Some(Self::Delta),
            "sam" => Some(Self::Sam),
            _ => None,
        }
    }
    pub fn validate_exact_output(self) -> Result<(), String> {
        match self {
            Self::Delta | Self::Sam => Err("SAM/delta output is not implemented correctly for this exact-match engine; use matches or paf".into()),
            _ => Ok(()),
        }
    }
}

/// Write record-local coordinates and preserve strand identity. PAF MAPQ 255
/// means unavailable; no confidence estimate is inferred from an exact seed.
pub fn write_record_matches<W: Write>(
    out: &mut W,
    matches: &[RecordMatch],
    references: &[FastaRecord],
    query: &FastaRecord,
    format: OutputFormat,
) -> io::Result<()> {
    format
        .validate_exact_output()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    if format == OutputFormat::Default {
        writeln!(out, "> {}", query.id)?;
    }
    for m in matches {
        let reference = references.get(m.reference_index).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "Invalid reference index")
        })?;
        if m.len == 0
            || m.ref_pos
                .checked_add(m.len)
                .is_none_or(|end| end > reference.sequence.len())
            || m.query_pos
                .checked_add(m.len)
                .is_none_or(|end| end > query.sequence.len())
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Match outside record boundaries",
            ));
        }
        match format {
            OutputFormat::Default => writeln!(
                out,
                "  Ref: {}:{}  Query: {}  Len: {}  Strand: {}",
                reference.id,
                m.ref_pos + 1,
                m.query_pos + 1,
                m.len,
                m.strand.symbol()
            )?,
            OutputFormat::Matches => writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}",
                reference.id,
                query.id,
                m.strand.symbol(),
                m.ref_pos,
                m.query_pos,
                m.len
            )?,
            OutputFormat::Paf => writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255",
                query.id,
                query.sequence.len(),
                m.query_pos,
                m.query_pos + m.len,
                m.strand.symbol(),
                reference.id,
                reference.sequence.len(),
                m.ref_pos,
                m.ref_pos + m.len,
                m.len,
                m.len
            )?,
            OutputFormat::Delta | OutputFormat::Sam => unreachable!("validated above"),
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Strand;
    #[test]
    fn paf_reverse_has_local_coordinates_strand_and_unknown_mapq() {
        let references = vec![FastaRecord {
            id: "r".into(),
            sequence: b"AACCGGTA".to_vec(),
        }];
        let query = FastaRecord {
            id: "q".into(),
            sequence: b"CCTACCGGTTG".to_vec(),
        };
        let matches = vec![RecordMatch {
            reference_index: 0,
            ref_pos: 0,
            query_pos: 2,
            len: 8,
            strand: Strand::Reverse,
        }];
        let mut out = Vec::new();
        write_record_matches(&mut out, &matches, &references, &query, OutputFormat::Paf).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "q\t11\t2\t10\t-\tr\t8\t0\t8\t8\t8\t255\n"
        );
    }
}
