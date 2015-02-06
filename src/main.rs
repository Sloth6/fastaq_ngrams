#![feature(collections,path,io,core)]
use std::old_io::File;
use std::os;
use std::old_io::BufferedReader;
use std::collections::HashMap;

#[derive(Hash)]
#[derive(PartialEq)]
#[derive(Eq)]
#[derive(Clone,Debug)]

enum Nucleotide { A, C, T, G, }

impl Nucleotide {
  fn stringify(&self) -> String {
    match *self {
      Nucleotide::A => "A".to_string(),
      Nucleotide::T => 'T'.to_string(),
      Nucleotide::C => 'C'.to_string(),
      Nucleotide::G => 'G'.to_string(),
    }
  }
}

type DNA = Vec<Nucleotide>;

impl DNA {
  fn stringify(&self) -> String {
    let mut s = "<".to_string();
    for n in self.iter() {
      s.push_str(&*n.stringify())
    }
    s.push_str(">");
    s
  }
}

struct Record {
  name: String,
  seq: Vec<Nucleotide>,
}

impl Record {
  fn stringify(&self) -> String {
    let mut s = self.name.to_string();
    s.push_str(": <");
    for n in self.seq.iter() {
      s.push_str(&*n.stringify());
    }
    s.push_str(">");
    s
  }
}

fn dna_to_String(DNA: &Vec<Nucleotide>) -> String {
  let mut s = "<".to_string();
  for n in DNA.iter() {
    s.push_str(&*n.stringify());
  }
  s.push_str(">");
  s
}

fn load_fastq(path: String) -> Vec<Record> {
  let path = Path::new(path);

  // Open the path in read-only mode, returns `IoResult<File>`
  let file = match File::open(&path) {
    // The `desc` field of `IoError` is a string that describes the error
    Err(why) => panic!("couldn't open {}: {}", path.display(), why.desc),
    Ok(file) => file,
  };

  let mut reader = BufferedReader::new(file);
  let mut lines = reader.lines().filter_map(|result| result.ok());

  let mut sequences: Vec<Record> = vec![];
  // let mut count: usize = 0;
  loop {
    let l_a = lines.next();
    let l_b = lines.next();
    let l_c = lines.next();
    let l_d = lines.next();
    match ( l_a ,l_b ,l_c ,l_d ) {
      (Some(name), Some(seq), _, _) => {
        let mut nuc_seq: Vec<Nucleotide> = vec![];
        
        for c in seq.chars() {
          match c {
            'A' => nuc_seq.push(Nucleotide::A),
            'T' => nuc_seq.push(Nucleotide::T),
            'C' => nuc_seq.push(Nucleotide::C),
            'G' => nuc_seq.push(Nucleotide::G),
            _ => (),
          }
        }
        sequences.push(Record { name: name, seq: nuc_seq });
        if sequences.len() >= 1000000 {
          break;
        }
      },
      _ => { break }
    }
  }
  println!("Loaded File");
  sequences
}

fn main() {
  
  // let sequences = load_fastq("./src/test.fastq".to_string());
  let sequences = load_fastq("/Users/joelsimon/Documents/Immufind/Seq Clustering/data/SRR015423_1.filt.fastq".to_string());
  let k: usize = 20;
  let mut n_grams: HashMap<Vec<Nucleotide>, usize> = HashMap::new();
  let mut max: usize = 0;
  for record in sequences.iter() {
    for i in range(0, (*record).seq.len() - k +1) {
      let mut subseq: Vec<Nucleotide> = vec![];
      for j in range(i, i+k) {
        subseq.push((*record).seq[j].clone());
      }
      let mut fuck = false;
      match n_grams.get_mut(&subseq) {
        Some(n) => {
          *n += 1;
          if *n > max {
            max = *n;
          }
        },
        None => {
          fuck = true;
        },
      }
      if fuck == true {
        n_grams.insert(subseq, 1);
      } 
    }
  }
  println!("{:?}", max);
  for (k, v) in n_grams.drain() {
    if v > 2000 {
      println!("{:?}{:?}", dna_to_String(&k), v );  
    }
  }  
}
      // println!("  {:?}", subRecord.stringify());

// fn levenshtein_distance(word1: &str, word2: &str) -> uint {
//   let word1_length = word1.len();
//   let word2_length = word2.len();
 
//   let mut matrix : Vec<Vec<uint>> = vec![vec![]];
 
//   for i in range(1, word1_length+1) { matrix[0].push(i); }
//   for j in range(1, word2_length+1) { matrix.push(vec![j]); }


//   for j in range(1, word2_length) {
//     for i in range(1, word1_length) {
//       let x: uint = if word1.char_at(i) == word2.char_at(j) {
//         matrix[j-1][i-1]
//       }
//       else {
//         let min_distance = [matrix[j][i-1], matrix[j-1][i], matrix[j-1][i-1]];
//         *min_distance.iter().min().unwrap() + 1
//       };
 
//       matrix[j].push(x);
//     }
//   }
 
//   matrix[word2_length-1][word1_length-1]
// }

// fn distance_matrix(seqs: Vec<Record>) -> Vec<Vec<uint>> {
//  let mut matrix : Vec<Vec<uint>> = vec![vec![]];
//  let n = seqs.len();

//  for _ in range(1, n) {
//    matrix.push(vec![]);
//  }
  
//  for i in range(0, n) {
//    for j in range(0, i+1) {
//      let d = levenshtein_distance(&seqs[i].seq, &seqs[j].seq);
//      matrix[i].push(d);
//    }
//  }

//  return matrix;
// }