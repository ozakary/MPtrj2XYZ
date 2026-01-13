#!/usr/bin/env python3
import ijson
from pathlib import Path
from typing import Set, Dict, Tuple
import argparse
from tqdm import tqdm
from ase import Atoms
from ase.io import write
import numpy as np
import time


class MPtrjExtractor:
    
    def __init__(self, target_elements: Set[str] = None, output_dir: str = "output", strategies: list = None):
        self.target_elements = target_elements
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.convert_all = (target_elements is None)
        
        if strategies is None:
            self.strategies = ['any', 'only', 'all']
        else:
            self.strategies = strategies
        
        self.stats = {strategy: {'count': 0, 'file': None} for strategy in self.strategies}
        self.stats['total_processed'] = 0
        self.stats['start_time'] = time.time()
        
        self._open_output_files()
    
    def _open_output_files(self):
        for strategy in self.strategies:
            filepath = self.output_dir / f"structures_{strategy}.xyz"
            self.stats[strategy]['file'] = open(filepath, 'w')
    
    def _close_output_files(self):
        for strategy in self.strategies:
            if self.stats[strategy]['file']:
                self.stats[strategy]['file'].close()
    
    def _parse_structure(self, structure_data: Dict) -> Tuple[Atoms, Set[str]]:
        structure = structure_data['structure']
        lattice_matrix = structure['lattice']['matrix']
        cell = np.array(lattice_matrix)
        
        symbols = []
        positions = []
        elements_in_structure = set()
        
        for site in structure['sites']:
            element = site['species'][0]['element']
            symbols.append(element)
            elements_in_structure.add(element)
            positions.append(site['xyz'])
        
        atoms = Atoms(symbols=symbols, positions=np.array(positions), cell=cell, pbc=True)
        return atoms, elements_in_structure
    
    def _check_filters(self, elements: Set[str]) -> Dict[str, bool]:
        if self.convert_all:
            # Convert all structures - always pass
            return {strategy: True for strategy in self.strategies}
        
        return {
            'any': len(elements & self.target_elements) > 0,
            'only': elements.issubset(self.target_elements) and len(elements) > 0,
            'all': self.target_elements.issubset(elements)
        }
    
    def _write_structure_to_xyz(self, atoms: Atoms, structure_id: str, strategy: str, 
                                structure_data: Dict, metadata: Dict = None):
        """Write structure to extended XYZ file with forces, energy, and stress."""
        
        # Get forces
        forces = structure_data.get('force', None)
        if forces is not None:
            forces_array = np.array(forces)
            atoms.set_array('forces', forces_array)
        
        # Get stress tensor (6 components or 3x3 matrix)
        stress = structure_data.get('stress', None)
        if stress is not None:
            stress_array = np.array(stress, dtype=float)  # Ensure float type to handle Decimal
            if stress_array.shape == (3, 3):
                # Already 3x3, just flatten to 9 components
                stress_full = stress_array.flatten().tolist()
            else:
                # Convert from Voigt notation (6 components: xx yy zz yz xz xy)
                # to full 3x3 matrix (9 components: xx xy xz yx yy yz zx zy zz)
                stress_voigt = stress_array.flatten()
                if len(stress_voigt) == 6:
                    sxx, syy, szz, syz, sxz, sxy = stress_voigt
                    # Build symmetric 3x3 matrix and flatten
                    stress_full = [
                        sxx, sxy, sxz,
                        sxy, syy, syz,
                        sxz, syz, szz
                    ]
                else:
                    stress_full = stress_voigt.tolist()
            
            atoms.info['stress'] = ' '.join(map(str, stress_full))
        
        # Get energies - ensure they are float scalars
        corrected_energy = structure_data.get('corrected_total_energy', None)
        uncorrected_energy = structure_data.get('uncorrected_total_energy', None)
        
        if corrected_energy is not None:
            atoms.info['free_energy'] = float(corrected_energy)
        if uncorrected_energy is not None:
            atoms.info['energy'] = float(uncorrected_energy)
        
        # Write to file using ASE extended XYZ format
        write(self.stats[strategy]['file'], atoms, format='extxyz')
        self.stats[strategy]['count'] += 1
    
    def process_json_file(self, json_filepath: str, chunk_size: int = 1000):
        print(f"Processing: {json_filepath}")
        if self.convert_all:
            print(f"Mode: Converting ALL structures (no element filtering)")
        else:
            print(f"Target elements: {sorted(self.target_elements)}")
        print(f"Output directory: {self.output_dir}")
        print(f"Strategies: {', '.join(self.strategies)}")
        
        if not self.convert_all:
            print("\nStrategy definitions:")
            print("  'any'  : At least one target element (can have others)")
            print("  'only' : Only target elements (no other elements)")
            print("  'all'  : All target elements must be present")
        print()
        
        with open(json_filepath, 'rb') as f:
            parser = ijson.kvitems(f, '')
            
            with tqdm(desc="Processing", unit=" struct", 
                     bar_format='{desc}: {n_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}') as pbar:
                
                for mp_id, mp_data in parser:
                    for structure_id, structure_data in mp_data.items():
                        try:
                            atoms, elements = self._parse_structure(structure_data)
                            filters_passed = self._check_filters(elements)
                            
                            metadata = {
                                'corrected_total_energy': structure_data.get('corrected_total_energy'),
                                'mp_id': structure_data.get('mp_id', mp_id)
                            }
                            
                            full_id = f"{mp_id}_{structure_id}"
                            for strategy in self.strategies:
                                if filters_passed[strategy]:
                                    self._write_structure_to_xyz(atoms, full_id, strategy, structure_data, metadata)
                            
                            self.stats['total_processed'] += 1
                            
                            if self.stats['total_processed'] % chunk_size == 0:
                                pbar.update(chunk_size)
                                
                                elapsed = time.time() - self.stats['start_time']
                                rate = self.stats['total_processed'] / elapsed
                                
                                postfix_dict = {'rate': f'{rate:.1f}/s'}
                                for strategy in self.strategies:
                                    postfix_dict[strategy] = self.stats[strategy]['count']
                                
                                pbar.set_postfix(postfix_dict)
                                
                        except Exception as e:
                            print(f"\nError: {mp_id}/{structure_id}: {e}")
                            continue
                
                pbar.update(self.stats['total_processed'] % chunk_size)
    
    def print_statistics(self):
        elapsed = time.time() - self.stats['start_time']
        
        print("\n" + "="*70)
        print("EXTRACTION STATISTICS")
        print("="*70)
        print(f"Total structures processed: {self.stats['total_processed']:,}")
        print(f"Total time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
        print(f"Processing rate: {self.stats['total_processed']/elapsed:.1f} structures/second")
        print()
        print("Extracted structures by strategy:")
        print("-" * 70)
        
        for strategy in self.strategies:
            count = self.stats[strategy]['count']
            percentage = (count / self.stats['total_processed'] * 100 
                         if self.stats['total_processed'] > 0 else 0)
            filepath = self.output_dir / f"structures_{strategy}.xyz"
            size_mb = filepath.stat().st_size / (1024 * 1024) if filepath.exists() else 0
            
            print(f"  {strategy:6s}: {count:>10,} structures ({percentage:>5.2f}%) | {size_mb:>8.2f} MB")
        
        print("="*70)
    
    def close(self):
        self._close_output_files()


def main():
    parser = argparse.ArgumentParser(
        description="Extract structures from MPtrj dataset based on element filtering"
    )
    
    parser.add_argument('json_file', type=str, help='Path to MPtrj JSON file')
    parser.add_argument('-s', '--strategy', nargs='+', choices=['any', 'only', 'all'],
                       default=['any', 'only', 'all'], help='Filter strategies (default: all)')
    parser.add_argument('-e', '--elements', nargs='+', default=['H', 'C', 'N', 'O', 'F', 'S', 'Cl', 'Xe'],
                       help='Target elements (default: H C N O F S Cl Xe). Use "all" to convert entire dataset')
    parser.add_argument('-o', '--output', type=str, default='mptrj_extracted',
                       help='Output directory (default: mptrj_extracted)')
    parser.add_argument('-c', '--chunk-size', type=int, default=1000,
                       help='Progress update frequency (default: 1000)')
    
    args = parser.parse_args()
    
    if not Path(args.json_file).exists():
        print(f"Error: File '{args.json_file}' not found!")
        return 1
    
    target_elements = set(args.elements)
    
    # Check if user wants all elements
    if 'all' in target_elements or args.elements == ['all']:
        # Use a special flag to bypass filtering
        target_elements = None
        strategies = ['all']  # Only use 'all' strategy for full conversion
    else:
        strategies = args.strategy
    
    extractor = MPtrjExtractor(target_elements, args.output, strategies)
    
    try:
        extractor.process_json_file(args.json_file, args.chunk_size)
        extractor.print_statistics()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user!")
        extractor.print_statistics()
        return 1
    except Exception as e:
        print(f"\nFatal error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    finally:
        extractor.close()
    
    print(f"\nOutput files saved to: {extractor.output_dir.absolute()}")
    return 0


if __name__ == "__main__":
    exit(main())
