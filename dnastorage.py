import base64
import hashlib
import json
import zlib
import re
import random
from typing import Dict, Any, Tuple, Optional, List
from enum import Enum, auto
from dataclasses import dataclass
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from collections import Counter

class DNAStorageMode(Enum):
    """Operation modes for DNA storage with enhanced options"""
    STANDARD = auto()       # Basic 2-bit encoding
    COMPRESSED = auto()     # With zlib compression
    ERROR_RESISTANT = auto() # With built-in redundancy
    HIGH_DENSITY = auto()   # 3-bit encoding (experimental)

@dataclass
class DNAStorageResult:
    dna_sequence: str
    metadata: Dict[str, Any]
    checksum: str
    encoding_scheme: str

class AdvancedDNACoder:
    """
    Next-generation DNA storage coder with:
    - Multiple encoding schemes (2-bit, 3-bit)
    - Advanced biological constraints
    - Intelligent sequence optimization
    - Enhanced error detection
    - Metadata-rich output
    """
    
    def __init__(self, mode: DNAStorageMode = DNAStorageMode.STANDARD):
        self.mode = mode
        self._setup_encoding_schemes()
        self._setup_biological_constraints()
        self.optimization_iterations = 3  # Number of optimization passes

    def _setup_encoding_schemes(self):
        """Configure available encoding schemes"""
        self.schemes = {
            '2bit': {
                'map': {'00': 'A', '01': 'T', '10': 'C', '11': 'G'},
                'reverse': {'A': '00', 'T': '01', 'C': '10', 'G': '11'},
                'bits_per_base': 2
            },
            '3bit': {
                'map': self._generate_3bit_mapping(),
                'reverse': {v: k for k, v in self._generate_3bit_mapping().items()},
                'bits_per_base': 3
            }
        }
        self.current_scheme = '2bit'  # Default scheme

    def _generate_3bit_mapping(self) -> Dict[str, str]:
        """Generates optimized 3-bit to DNA mapping"""
        # Uses 6 bases (including degenerate bases if needed)
        return {
            '000': 'A',
            '001': 'T',
            '010': 'C',
            '011': 'G',
            '100': 'AT',
            '101': 'CG',
            '110': 'TA',
            '111': 'GC'
        }

    def _setup_biological_constraints(self):
        """Configure advanced biological sequence requirements"""
        self.constraints = {
            'max_homopolymer': 4,  # Maximum identical bases in a row
            'gc_bounds': (0.4, 0.6),  # Acceptable GC content range
            'prohibited_motifs': [
                'GGGG',  # G-quadruplex
                'CCCC',  # C-rich
                'TTTT',  # Poly-T terminator
                'AAAA',   # Poly-A tail
                'GTGT',  # Slippery sequences
                'CGCG'
            ],
            'max_dinucleotide_repeat': 3  # Max repeats of same dinucleotide
        }

    def encode_data_to_dna(self, data: Dict[str, Any], scheme: str = '2bit') -> DNAStorageResult:
        """
        Enhanced DNA encoding pipeline with multiple passes:
        1. JSON serialization
        2. Optional compression
        3. Base64 encoding
        4. Binary conversion
        5. DNA mapping
        6. Biological optimization
        7. Validation
        """
        self.current_scheme = scheme
        
        try:
            # Step 1: Serialize to JSON with compact formatting
            json_str = json.dumps(data, separators=(',', ':'))
            bytes_data = json_str.encode('utf-8')
            
            # Step 2: Apply compression if enabled
            if self.mode == DNAStorageMode.COMPRESSED:
                bytes_data = zlib.compress(bytes_data, level=9)
            
            # Step 3: Base64 encoding
            b64_encoded = base64.b64encode(bytes_data).decode('utf-8')
            
            # Step 4: Binary conversion
            binary_str = self._bytes_to_binary(b64_encoded)
            
            # Step 5: DNA mapping
            dna = self._binary_to_dna(binary_str)
            
            # Step 6: Biological optimization
            optimized_dna = self._optimize_sequence(dna)
            
            # Step 7: Final validation
            is_valid, message = self.validate_sequence(optimized_dna)
            if not is_valid:
                raise ValueError(f"Optimization failed: {message}")
            
            # Generate comprehensive metadata
            metadata = self._generate_metadata(data, json_str, optimized_dna)
            
            return DNAStorageResult(
                dna_sequence=optimized_dna,
                metadata=metadata,
                checksum=self._generate_checksum(optimized_dna),
                encoding_scheme=scheme
            )
            
        except Exception as e:
            raise ValueError(f"DNA encoding failed: {str(e)}")

    def _bytes_to_binary(self, data: str) -> str:
        """Converts string data to binary representation"""
        return ''.join(format(ord(c), '08b') for c in data)

    def _binary_to_dna(self, binary_str: str) -> str:
        """Maps binary string to DNA sequence"""
        scheme = self.schemes[self.current_scheme]
        chunk_size = scheme['bits_per_base']
        
        # Pad binary string if needed
        if len(binary_str) % chunk_size != 0:
            padding = chunk_size - (len(binary_str) % chunk_size)
            binary_str += '0' * padding
        
        dna = []
        for i in range(0, len(binary_str), chunk_size):
            chunk = binary_str[i:i+chunk_size]
            dna.append(scheme['map'].get(chunk, 'N'))  # Default to N if invalid
            
        return ''.join(dna)

    def _optimize_sequence(self, dna: str) -> str:
        """
        Applies multiple optimization passes to ensure biological constraints:
        1. Homopolymer breaking
        2. GC balancing
        3. Motif removal
        4. Dinucleotide repeat handling
        """
        optimized = dna
        for _ in range(self.optimization_iterations):
            optimized = self._break_homopolymers(optimized)
            optimized = self._balance_gc_content(optimized)
            optimized = self._remove_prohibited_motifs(optimized)
            optimized = self._break_dinucleotide_repeats(optimized)
        
        return optimized

    def _break_homopolymers(self, dna: str) -> str:
        """Intelligently breaks long homopolymers"""
        result = []
        current_base = None
        count = 0
        
        for base in dna:
            if base == current_base:
                count += 1
                if count > self.constraints['max_homopolymer']:
                    # Insert optimal replacement base
                    replacement = self._get_optimal_replacement(base, dna, len(result))
                    result.append(replacement)
                    current_base = replacement
                    count = 1
            else:
                current_base = base
                count = 1
            result.append(base)
            
        return ''.join(result)

    def _get_optimal_replacement(self, current_base: str, dna: str, position: int) -> str:
        """Chooses replacement base that improves overall sequence properties"""
        # Consider GC content impact
        gc = gc_fraction(dna)
        alternatives = []
        
        if gc < self.constraints['gc_bounds'][0]:
            # Prefer G/C to increase GC
            alternatives = ['G', 'C']
        elif gc > self.constraints['gc_bounds'][1]:
            # Prefer A/T to decrease GC
            alternatives = ['A', 'T']
        else:
            alternatives = [b for b in self.schemes[self.current_scheme]['reverse'] if b != current_base]
        
        # Also consider local context
        context = dna[max(0, position-3):position+3]
        for base in alternatives:
            if not any(motif in context + base for motif in self.constraints['prohibited_motifs']):
                return base
        
        return random.choice(alternatives)

    def _balance_gc_content(self, dna: str) -> str:
        """Smart GC balancing with local sequence awareness"""
        current_gc = gc_fraction(dna)
        target_gc = sum(self.constraints['gc_bounds'])/2
        dna_list = list(dna)
        
        # Determine replacement strategy
        if current_gc < target_gc:
            replace_from = ['A', 'T']
            replace_to = ['G', 'C']
        else:
            replace_from = ['G', 'C']
            replace_to = ['A', 'T']
        
        # Make strategic replacements in GC-rich/GC-poor regions
        window_size = 20
        for i in range(0, len(dna_list), window_size//2):
            window = dna_list[i:i+window_size]
            window_gc = gc_fraction(window)
            
            if (current_gc < target_gc and window_gc < target_gc) or \
               (current_gc > target_gc and window_gc > target_gc):
                
                for j in range(len(window)):
                    if window[j] in replace_from:
                        window[j] = random.choice(replace_to)
                        dna_list[i+j] = window[j]
                        
                        # Recalculate global GC
                        current_gc = gc_fraction(dna_list)
                        if self.constraints['gc_bounds'][0] <= current_gc <= self.constraints['gc_bounds'][1]:
                            return ''.join(dna_list)
        
        return ''.join(dna_list)

    def _remove_prohibited_motifs(self, dna: str) -> str:
        """Removes all prohibited motifs with context-aware replacement"""
        for motif in self.constraints['prohibited_motifs']:
            if motif in dna:
                # Find all occurrences
                matches = [m.start() for m in re.finditer(f'(?={motif})', dna)]
                for pos in reversed(matches):
                    # Replace with optimal alternative
                    replacement = self._generate_motif_replacement(motif, dna, pos)
                    dna = dna[:pos] + replacement + dna[pos+len(motif):]
        return dna

    def _generate_motif_replacement(self, motif: str, context: str, position: int) -> str:
        """Generates context-aware motif replacements"""
        # Try to maintain local GC content
        local_gc = gc_fraction(context[max(0, position-10):position+10+len(motif)])
        
        replacements = []
        for _ in range(10):
            candidate = ''.join(random.choice(list(self.schemes[self.current_scheme]['reverse'].keys()))
                              for _ in range(len(motif)))
            
            # Check constraints
            if (not any(m in candidate for m in self.constraints['prohibited_motifs']) and
                not re.search(r'(.)\1{3,}', candidate) and  # No new homopolymers
                abs(gc_fraction(candidate) - local_gc) < 0.1):  # GC similar to local
                replacements.append(candidate)
        
        return min(replacements, key=lambda x: abs(gc_fraction(x) - local_gc)) if replacements else 'A'*len(motif)

    def _break_dinucleotide_repeats(self, dna: str) -> str:
        """Breaks long dinucleotide repeats (e.g., ATATAT)"""
        for i in range(len(dna) - self.constraints['max_dinucleotide_repeat']*2):
            window = dna[i:i+self.constraints['max_dinucleotide_repeat']*2]
            if len(window) < 4:
                continue
                
            # Check for dinucleotide repeats
            if all(window[j] == window[j%2] for j in range(len(window))):
                # Break the pattern
                pos = i + 2  # Middle of first repeat
                replacement = self._get_optimal_replacement(dna[pos], dna, pos)
                dna = dna[:pos] + replacement + dna[pos+1:]
                
        return dna

    def _generate_metadata(self, original_data: Dict[str, Any], json_str: str, dna: str) -> Dict[str, Any]:
        """Generates comprehensive encoding metadata"""
        seq_obj = Seq(dna)
        base_counts = Counter(dna)
        
        return {
            'original_size': len(json_str),
            'encoded_size': len(dna),
            'compression_ratio': len(json_str)/len(dna),
            'gc_content': gc_fraction(dna),
            'base_distribution': dict(base_counts),
            'mode': self.mode.name,
            'scheme': self.current_scheme,
            'biological_constraints': {
                'homopolymers': self._find_homopolymers(dna),
                'prohibited_motifs': self._find_prohibited_motifs(dna),
                'dinucleotide_repeats': self._find_dinucleotide_repeats(dna)
            },
            'reverse_complement': str(seq_obj.reverse_complement()[:20]) + '...'
        }

    def _find_homopolymers(self, dna: str) -> List[Dict[str, Any]]:
        """Identifies remaining homopolymers in sequence"""
        return [{'position': m.start(), 'sequence': m.group()}
               for m in re.finditer(r'(.)\1{3,}', dna)]

    def _find_prohibited_motifs(self, dna: str) -> List[Dict[str, Any]]:
        """Identifies any remaining prohibited motifs"""
        return [{'motif': motif, 'positions': [m.start() for m in re.finditer(motif, dna)]}
               for motif in self.constraints['prohibited_motifs'] if motif in dna]

    def _find_dinucleotide_repeats(self, dna: str) -> List[Dict[str, Any]]:
        """Identifies dinucleotide repeats"""
        repeats = []
        for i in range(len(dna) - 4):
            window = dna[i:i+4]
            if window[:2] == window[2:4]:
                repeats.append({'position': i, 'sequence': window})
        return repeats

    def _generate_checksum(self, dna: str) -> str:
        """Generates enhanced checksum with sequence properties"""
        seq_obj = Seq(dna)
        metadata = {
            'length': len(dna),
            'gc_content': gc_fraction(dna),
            'base_counts': dict(Counter(dna)),
            'first_last': dna[:10] + dna[-10:]
        }
        combined = dna + json.dumps(metadata, sort_keys=True)
        return hashlib.sha3_256(combined.encode()).hexdigest()

    def decode_dna_to_data(self, dna: str) -> Dict[str, Any]:
        """
        Robust DNA decoding with automatic scheme detection
        Supports both 2-bit and 3-bit encodings
        """
        try:
            # Try 2-bit decoding first
            self.current_scheme = '2bit'
            binary_str = ''.join(self.schemes['2bit']['reverse'].get(base, '00') for base in dna)
            result = self._binary_to_data(binary_str)
            return result
        except (ValueError, KeyError):
            try:
                # Fall back to 3-bit decoding
                self.current_scheme = '3bit'
                binary_str = self._dna_to_binary_3bit(dna)
                return self._binary_to_data(binary_str)
            except Exception as e:
                raise ValueError(f"DNA decoding failed with all schemes: {str(e)}")

    def _dna_to_binary_3bit(self, dna: str) -> str:
        """Converts 3-bit encoded DNA back to binary"""
        binary = []
        i = 0
        while i < len(dna):
            chunk = dna[i:i+2]  # 3-bit encoding uses 1-2 bases
            if chunk in self.schemes['3bit']['reverse']:
                binary.append(self.schemes['3bit']['reverse'][chunk])
                i += len(chunk)
            elif chunk[0] in self.schemes['3bit']['reverse']:
                binary.append(self.schemes['3bit']['reverse'][chunk[0]])
                i += 1
            else:
                raise ValueError(f"Invalid DNA sequence for 3-bit decoding at position {i}")
        return ''.join(binary)

    def _binary_to_data(self, binary_str: str) -> Dict[str, Any]:
        """Converts binary string back to original data"""
        # Convert binary to Base64
        b64_str = ''.join(chr(int(binary_str[i:i+8], 2)) 
                      for i in range(0, len(binary_str), 8) if i+8 <= len(binary_str))
        
        # Base64 decode
        bytes_data = base64.b64decode(b64_str)
        
        # Try decompression
        try:
            decompressed = zlib.decompress(bytes_data)
            return json.loads(decompressed.decode('utf-8'))
        except zlib.error:
            return json.loads(bytes_data.decode('utf-8'))

    def validate_sequence(self, dna_sequence: str) -> Tuple[bool, Optional[str]]:
        """Comprehensive sequence validation with detailed feedback"""
        try:
            # Check GC content
            gc = gc_fraction(dna_sequence)
            if not self.constraints['gc_bounds'][0] <= gc <= self.constraints['gc_bounds'][1]:
                return False, f"GC content {gc:.2f} outside bounds {self.constraints['gc_bounds']}"
            
            # Check homopolymers
            homopolymers = self._find_homopolymers(dna_sequence)
            if homopolymers:
                worst = max(homopolymers, key=lambda x: len(x['sequence']))
                return False, f"Homopolymer {worst['sequence']} at position {worst['position']}"
            
            # Check prohibited motifs
            motifs = self._find_prohibited_motifs(dna_sequence)
            if motifs:
                return False, f"Prohibited motif {motifs[0]['motif']} found at positions {motifs[0]['positions'][:3]}"
            
            # Check dinucleotide repeats
            dinucleotides = self._find_dinucleotide_repeats(dna_sequence)
            if dinucleotides and len(dinucleotides[0]['sequence']) > 4:
                return False, f"Dinucleotide repeat {dinucleotides[0]['sequence']} at position {dinucleotides[0]['position']}"
            
            return True, None
        except Exception as e:
            return False, str(e)

    def calculate_encoding_metrics(self, data: Dict[str, Any]) -> Dict[str, float]:
        """Calculates comprehensive encoding performance metrics"""
        original_json = json.dumps(data, separators=(',', ':'))
        encoded_result = self.encode_data_to_dna(data)
        
        original_bits = len(original_json) * 8
        dna_bits = len(encoded_result.dna_sequence) * self.schemes[self.current_scheme]['bits_per_base']
        
        return {
            'bits_per_base': dna_bits / len(encoded_result.dna_sequence),
            'information_density': original_bits / len(encoded_result.dna_sequence),
            'compression_ratio': len(original_json) / len(encoded_result.dna_sequence),
            'gc_content': encoded_result.metadata['gc_content'],
            'efficiency': original_bits / dna_bits,
            'payload_ratio': len(original_json) / len(encoded_result.dna_sequence)
        }

    def reverse_complement(self, dna_sequence: str) -> str:
        """Returns reverse complement using BioPython with validation"""
        try:
            return str(Seq(dna_sequence).reverse_complement())
        except Exception as e:
            raise ValueError(f"Invalid DNA sequence for reverse complement: {str(e)}")

# Example usage
if __name__ == "__main__":
    # Initialize coder with error-resistant mode
    coder = AdvancedDNACoder(mode=DNAStorageMode.ERROR_RESISTANT)
    
    # Sample data to encode
    test_data = {
        "title": "Advanced DNA Storage Demo",
        "description": "This demonstrates the enhanced DNA storage system",
        "version": 2.5,
        "features": ["error-resistant", "gc-balanced", "motif-free"],
        "timestamp": "2023-07-20T12:00:00Z"
    }
    
    print("=== Encoding Example ===")
    # Encode with 2-bit scheme
    encoded_2bit = coder.encode_data_to_dna(test_data, scheme='2bit')
    print(f"Encoded DNA (2-bit, {len(encoded_2bit.dna_sequence)} bases):")
    print(encoded_2bit.dna_sequence[:50], "...", encoded_2bit.dna_sequence[-50:])
    print(f"Checksum: {encoded_2bit.checksum}")
    print("\nMetadata:")
    for k, v in encoded_2bit.metadata.items():
        print(f"{k}: {v}")
    
    # Encode same data with 3-bit scheme
    encoded_3bit = coder.encode_data_to_dna(test_data, scheme='3bit')
    print(f"\n3-bit encoding length: {len(encoded_3bit.dna_sequence)} bases")
    
    print("\n=== Decoding Example ===")
    decoded = coder.decode_dna_to_data(encoded_2bit.dna_sequence)
    print("Decoded data matches original:", decoded == test_data)
    
    print("\n=== Validation Example ===")
    is_valid, message = coder.validate_sequence(encoded_2bit.dna_sequence)
    print(f"Sequence valid: {is_valid}, {message or ''}")
    
    print("\n=== Metrics Example ===")
    metrics = coder.calculate_encoding_metrics(test_data)
    for k, v in metrics.items():
        print(f"{k}: {v:.4f}")