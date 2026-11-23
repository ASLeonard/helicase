use super::*;
use fastx_gen::FastxGen;
use std::io::Cursor;

const CONFIG: Config = crate::ParserOptions::default().config();

fn count(mut parser: FastxParser<'_, CONFIG>) -> Result<(usize, usize)> {
    let (mut records, mut len) = (0, 0);
    while let Some(event) = parser.next() {
        event?;
        records += 1;
        len += parser.get_dna_string().len();
    }
    Ok((records, len))
}

fn count_stream(input: &[u8], threads: usize, chunk_len: usize) -> Result<(usize, usize)> {
    let mut parser =
        ParallelFastxParser::<CONFIG>::from_reader(Cursor::new(input)).threads(threads);
    parser.max_chunk_len = chunk_len;
    let counts = parser.for_each_with(
        || (0, 0),
        |acc, p| {
            let (records, len) = count(p)?;
            *acc = (acc.0 + records, acc.1 + len);
            Ok(())
        },
    )?;
    Ok(counts
        .into_iter()
        .fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1)))
}

#[test]
fn stream_chunks_match_sequential() {
    for gen_records in [FastxGen::fasta, FastxGen::fastq] {
        let records = |n, len| gen_records(n, len).with_seed(len as u64).generate_bytes();
        // Records much longer than the chunks, at the start, in the middle and at the end
        let lens = [
            (1, 20_000),
            (300, 50),
            (1, 50_000),
            (2, 30_000),
            (100, 10),
            (1, 40_000),
        ];
        let input: Vec<u8> = lens.iter().flat_map(|&(n, len)| records(n, len)).collect();
        let expected = count(FastxParser::from_slice(&input).unwrap()).unwrap();
        for threads in [1, 3] {
            for chunk_len in [100, 4096] {
                let res = count_stream(&input, threads, chunk_len).unwrap();
                assert_eq!(res, expected, "{threads} threads, chunks of {chunk_len}");
            }
        }
    }
    assert_eq!(count_stream(b"", 2, 100).unwrap(), (0, 0));
}

#[test]
fn stream_stops_on_error() {
    let input = FastxGen::fasta(1000, 100).generate_bytes();
    let mut parser = ParallelFastxParser::<CONFIG>::from_reader(Cursor::new(&input)).threads(2);
    parser.max_chunk_len = 100;
    let res = parser.for_each(|_| Err("".into()));
    assert!(res.is_err());
    assert!(count_stream(b"xyz\n".repeat(100).as_slice(), 2, 100).is_err());
}

#[test]
fn thread_file_follows_path() {
    let dir = std::env::temp_dir();
    let (a, b) = (
        dir.join("helicase_thread_file_a"),
        dir.join("helicase_thread_file_b"),
    );
    std::fs::write(&a, b"a").unwrap();
    std::fs::write(&b, b"b").unwrap();
    for (path, expected) in [(&a, b'a'), (&b, b'b'), (&a, b'a')] {
        let mut byte = [0];
        with_thread_file(path, |file| {
            file.seek(SeekFrom::Start(0))?;
            Ok(file.read_exact(&mut byte)?)
        })
        .unwrap();
        assert_eq!(byte[0], expected);
    }
    let _ = std::fs::remove_file(a);
    let _ = std::fs::remove_file(b);
}

#[test]
fn drop_cache_range_aligns_inward() {
    assert_eq!(aligned_drop_range(0, 0), None);
    assert_eq!(aligned_drop_range(128 << 10, 64 << 10), None);
    assert_eq!(aligned_drop_range(100, 64 << 10), None);
    assert_eq!(
        aligned_drop_range(0, 64 << 10),
        Some((0, 64 << 10))
    );
    assert_eq!(
        aligned_drop_range(1, (192 << 10) + 123),
        Some((64 << 10, 128 << 10))
    );
    assert_eq!(aligned_drop_range(u64::MAX, u64::MAX), None);
}

#[test]
fn drop_behind_file_reads_exact_bytes() {
    let path = std::env::temp_dir().join("helicase_drop_behind_file_test");
    let expected: Vec<u8> = (0..(200 << 10)).map(|i| (i * 13 + 7) as u8).collect();
    std::fs::write(&path, &expected).unwrap();
    let mut reader = DropBehindFile {
        file: File::open(&path).unwrap(),
        pos: 0,
        dropped_to: 0,
        enabled: true,
    };
    let mut actual = Vec::new();
    let mut chunk = [0u8; 4093];
    loop {
        let n = reader.read(&mut chunk).unwrap();
        if n == 0 {
            break;
        }
        actual.extend_from_slice(&chunk[..n]);
    }
    assert_eq!(actual, expected);
    let _ = std::fs::remove_file(path);
}
