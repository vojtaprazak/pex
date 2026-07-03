# -*- coding: utf-8 -*-
"""
iMOD_comfile.py

IMOD comfile parser/writer.

Design goals:
- Loss tolerant parsing of IMOD parameter formats
- Preserve original file text when unmodified (bytewise round-trip)
- Represent comfiles as ordered blocks to preserve structure and order
- Support safe mutation and deterministic rewriting
- Minimal dependencies, use os.path, compact code style
"""

from __future__ import print_function
import os

class IMOD_comfile:
    """
    Parse and manipulate IMOD comfiles.

    Usage:
        from iMOD_comfile import IMOD_comfile
        obj = IMOD_comfile(rec_dir, 'tilt.com', make_paths_absolute=False)
    """

    def __init__(self, rec_dir, com_name, make_paths_absolute=True, **kwargs):
        self.rec_dir = os.path.realpath(rec_dir)
        self.com_name = com_name
        self.make_paths_absolute = make_paths_absolute
        self.out_dir = kwargs.get('out_dir')
        # canonical linear representation
        self.blocks = []  # each block: {'header': str, 'params': {k:v}, 'separators': {k:sep}, 'footer': [str,...]}
        self.excludelist = []
        # bookkeeping
        self._original_text = ''
        self._parsed_dirty = False  # become True when user mutates data
        # read if file exists
        if os.path.isfile(os.path.join(self.rec_dir, self.com_name)):
            self.read_comfile()

    @property
    def dict(self):
        """Compatibility shim for older scripts expecting .dict access."""
        # Returns a combined dictionary of all parameters from all blocks
        combined = {}
        print('old IMOD_comfile api')
        for block in self.blocks:
            combined.update(block['params'])
        return combined

    # wrappers #####################
    def get_val(self, key):
        """
        Search all blocks for a key. 
        If not found, prints a block summary and raises KeyError.
        """
        for i, block in enumerate(self.blocks):
            if key in block['params']:
                return block['params'][key]
        
        # Error handling: Print summary before crashing
        print(f"\n[ERROR] Key '{key}' not found in {self.com_name}")
        print("-" * 40)
        self.pretty_print_blocks()
        print("-" * 40)
        
        raise KeyError(f"Parameter '{key}' is missing from {self.com_name}")

    def set_val(self, key, value, separator=None, auto_add=True):
        """
        Updates key if found (case-insensitive). 
        If not found and auto_add is True, adds to the first executable block.
        """
        target = key.upper()
        for i, block in enumerate(self.blocks):
            for k in list(block['params'].keys()):
                if k.upper() == target:
                    self.set_param(i, k, value, separator=separator)
                    return
        
        if auto_add:
            # Fallback: add to the block containing $tilt or $newstack (usually index 1)
            target_block = 1 if len(self.blocks) > 1 else 0
            self.set_param(target_block, key, value, separator=separator)
        else:
            raise KeyError(f"{key} not found and auto_add is False")

    def del_val(self, key):
        """
        Search all blocks for key and remove it.
        Marks object as dirty for reconstruction.
        """
        for i, block in enumerate(self.blocks):
            if key in block['params']:
                del block['params'][key]
                if key in block['separators']:
                    del block['separators'][key]
                self._parsed_dirty = True
                return

        print(f"\n[ERROR] Cannot delete '{key}'; key not found in {self.com_name}")
        self.pretty_print_blocks()
        raise KeyError(f"Key '{key}' not found in any block of {self.com_name}")

    def __contains__(self, key):
            """Allows usage: if 'ScaleShifts' in i_align:"""
            for block in self.blocks:
                if key in block['params']:
                    return True
            return False

    def ensure_val(self, key, value, block_idx=1, separator=','):
            """Updates the value if it exists, otherwise adds it to the specified block."""
            if key in self:
                self.set_val(key, value, separator=separator)
            else:
                # Fallback to adding it to the main program block (usually index 1)
                self.set_param(block_idx, key, value, separator=separator)
            
    # ---------------------------
    # low-level token parse helper
    # ---------------------------
    def parse_mixed_entry(self, val, path=False):
        """
        Parse IMOD-style parameter values.

        Returns (value, separator) where separator is ',' or ' ' or None.
        value is either scalar (int, float, str) or list of those, or '' for empty.
        """
        # Normalize input to token list
        if isinstance(val, (list, tuple)):
            tokens = []
            for v in val:
                # preserve existing tokens exactly if they're already separated
                s = str(v)
                # If original token contains commas, we keep them as one token for now
                if ',' in s and len(s.split()) == 1:
                    tokens.append(s)
                else:
                    tokens.extend(s.split())
        else:
            tokens = str(val).split()

        # detect top-level separator
        separator = None
        if len(tokens) == 1 and ',' in tokens[0]:
            tokens = tokens[0].split(',')
            separator = ','
        elif len(tokens) > 1:
            separator = ' '

        parsed = []
        for tok in tokens:
            if tok == '':
                continue
            # integer range of form a-b (only one dash allowed, negatives guarded)
            if tok.count('-') == 1:
                a, b = tok.split('-')
                # ensure both sides are ints (allow negative signs)
                if a.lstrip('-').isdigit() and b.lstrip('-').isdigit():
                    start = int(a)
                    end = int(b)
                    if start <= end:
                        for i in range(start, end + 1):
                            parsed.append(i)
                        continue
            # try int
            try:
                ival = int(tok)
                parsed.append(ival)
                continue
            except:
                pass
            # try float
            try:
                fval = float(tok)
                parsed.append(fval)
                continue
            except:
                pass
            # path handling only if requested (do not join absolute paths)
            if path and not os.path.isabs(tok):
                parsed.append(os.path.join(path, tok))
            else:
                parsed.append(tok)

        if len(parsed) == 0:
            return '', separator
        if len(parsed) == 1:
            return parsed[0], separator
        return parsed, separator

    # ---------------------------
    # parsing and representation
    # ---------------------------
    def read_comfile(self, com_name=False, path=False):
        """
        Parse comfile into self.blocks. Also save original text for exact round-trip if unchanged.
        """
        if com_name:
            self.com_name = com_name
        full = os.path.join(self.rec_dir, self.com_name)
        if not os.path.isfile(full):
            raise Exception("File not found. %s" % os.path.realpath(full))

        with open(full, 'r') as f:
            text = f.read()
        self._original_text = text
        lines = text.splitlines()

        # parse into blocks
        self.blocks = []
        current = None
        use_path = self.rec_dir if self.make_paths_absolute else False

        for raw in lines:
            line = raw.rstrip('\n')
            if line.strip() == '' or line.strip().startswith('#'):
                # preserve comments and blank lines inside a footer of current block if present
                if current is None:
                    # leading comment lines are represented as a block with header None
                    # if needed one could store global header; for now we preserve as a pseudo-block header
                    # we prefer to preserve exact original text by using _original_text when unchanged
                    continue
                else:
                    current.setdefault('footer', []).append(line)
                    continue

            if line.startswith('$'):
                # start new block
                if current is not None:
                    self.blocks.append(current)
                current = {'header': line, 'params': {}, 'separators': {}, 'footer': []}
                continue

            # ignore lines outside of blocks
            if current is None:
                continue

            # s = line.split()
            # key = s[0]
            # rest = s[1:] if len(s) > 1 else []
            # val, sep = self.parse_mixed_entry(rest, path=use_path)
            # current['params'][key] = val
            # current['separators'][key] = sep

            s = line.split()
            if not s: continue
            
            if len(s) == 1:
                current.setdefault('positional', []).append(s[0])
                continue

            key = s[0]
            rest = s[1:]
            val, sep = self.parse_mixed_entry(rest, path=use_path)
            current['params'][key] = val
            current['separators'][key] = sep

        
        if current is not None:
            self.blocks.append(current)

        # recompute excludelist (simple aggregation)
        self._merge_excludelists()
        self._parsed_dirty = False

    def _merge_excludelists(self):
        """Consolidates ExcludeList, EXCLUDELIST, etc., into self.excludelist and clears old keys."""
        keys_to_check = ['ExcludeList', 'EXCLUDELIST', 'ExcludeSections', 'EXCLUDELIST2']
        found = []
        for block in self.blocks:
            for k in keys_to_check:
                if k in block['params']:
                    v = block['params'][k]
                    if isinstance(v, (list, tuple)):
                        found.extend([int(x) for x in v])
                    else:
                        try:
                            found.append(int(v))
                        except (ValueError, TypeError):
                            pass
                    # IMPORTANT: Remove the key so it doesn't persist during write_comfile
                    del block['params'][k]
                    if k in block['separators']:
                        del block['separators'][k]
                    self._parsed_dirty = True
        
        self.excludelist = sorted(list(set(found)))
        return self.excludelist

    # ---------------------------
    # mutation helpers
    # ---------------------------
    def set_param(self, block_idx, key, value, separator=None):
        """
        Set parameter key in block number block_idx.
        Marks object as modified so write will reconstruct.
        """
        if block_idx < 0 or block_idx >= len(self.blocks):
            raise IndexError('block index out of range')
        self.blocks[block_idx]['params'][key] = value
        if separator is not None:
            self.blocks[block_idx]['separators'][key] = separator
        else:
            # if not specified, infer comma for lists
            if isinstance(value, (list, tuple)):
                self.blocks[block_idx]['separators'][key] = ','
            else:
                self.blocks[block_idx]['separators'][key] = ''
        self._parsed_dirty = True

    # ---------------------------
    # writing
    # ---------------------------
    def _val2str(self, val):
        """Return list of strings for a value."""
        if isinstance(val, (list, tuple)):
            out = []
            for v in val:
                if isinstance(v, float):
                    out.append(str(round(v, 3)))
                else:
                    out.append(str(v))
            return out
        if val is None:
            return []
        if isinstance(val, float):
            return [str(round(val, 3))]
        return [str(val)]

    def write_comfile(self, out_dir=None, change_name=False):
        """
        Write the comfile to disk.

        Behaviour depends on whether the comfile has been modified in memory:

        1) If the comfile has NOT been modified since it was read
           (i.e. no calls to set_param or other mutations),
           the original file contents are written back verbatim.
           This preserves bytewise identity, including comments,
           shell directives ($setenv, $if, etc.), spacing, and ordering.

        2) If the comfile HAS been modified,
           the file is reconstructed from the internal block representation.
           In this case:
             - All executable blocks ($ commands) are written in order
             - Parameters are written in a normalized form
             - A default '# IMOD command file' header is written
             - Exact original formatting is not guaranteed,
               but semantic equivalence is preserved

        Parameters
        ----------
        out_dir : str or None
            Directory to write the comfile into.
            If None, self.out_dir must be set.
        change_name : str or False
            Output filename. If False, the original comfile name is used.
            Use with care: IMOD expects specific comfile names
            (e.g. newst.com, tilt.com, align.com) in reconstruction workflows.

        Notes
        -----
        - This method does NOT execute the comfile.
        - This method does NOT attempt to interpret shell directives.
        - The caller is responsible for selecting which block(s) to execute.
        - For IMOD execution, prefer generating commands from specific blocks
          rather than relying on the first block in the file.

        Intended use cases
        ------------------
        - Copying comfiles without modification
        - Safely editing parameters (e.g. binning, output paths)
        - Generating reproducible, inspectable IMOD workflows

        See also
        --------
        set_param
        get_command_from_block
        pretty_print_blocks
        """
        if out_dir is None:
            out_dir = self.out_dir
        if out_dir is None:
            raise Exception('Set IMOD_comfile.out_dir or specify output directory.')
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        name = change_name if change_name else self.com_name
        out_path = os.path.join(out_dir, name)
        # if not mutated, preserve original
        if not self._parsed_dirty and self._original_text:
            with open(out_path, 'w') as f:
                f.write(self._original_text)
            return

        # else reconstruct
        with open(out_path, 'w') as f:
            # write a default header if none present originally
            f.write('# IMOD command file\n')

            for block in self.blocks:
                f.write(block['header'] + '\n')
                # Write positional arguments (0, filenames, etc.)
                for arg in block.get('positional', []):
                    f.write(arg + '\n')
                for key in block['params']:                    
            # for block in self.blocks:
            #     f.write(block['header'] + '\n')
            #     for key in block['params']:
                    sep = block['separators'].get(key)
                    if sep is None:
                        # default: comma for lists, empty for scalars
                        if isinstance(block['params'][key], (list, tuple)):
                            sep = ','
                        else:
                            sep = ''
                    vals = self._val2str(block['params'][key])
                    f.write('%s\t%s\n' % (key, sep.join(vals)))
                # footer lines
                for foot in block.get('footer', []):
                    f.write(foot + '\n')

    # ---------------------------
    # command construction
    # ---------------------------
    def get_command_list(self, exclude_keys=None, append_to_exclude_keys=None):
        """
        Build command list for subprocess from the first block.
        exclude_keys is list of keys to exclude.
        """
        if exclude_keys is None:
            exclude_keys = ['RADIAL', 'FalloffIsTrueSigma', 'ActionIfGPUFails', 'UseGPU', 'FakeSIRTiterations', 'SizeToOutputInXandY']
        if append_to_exclude_keys:
            exclude_keys = list(exclude_keys) + list(append_to_exclude_keys)

        if not self.blocks:
            return []

        head = self.blocks[0]['header']
        program = head.split()[0].lstrip('$')
        cmd = [program]
        for key in self.blocks[0]['params']:
            if key in exclude_keys:
                continue
            sep = self.blocks[0]['separators'].get(key)
            if sep is None:
                if isinstance(self.blocks[0]['params'][key], (list, tuple)):
                    sep = ','
                else:
                    sep = ''
            val = sep.join(self._val2str(self.blocks[0]['params'][key]))
            cmd.append('-%s' % key)
            cmd.append(val)
        return cmd

    def get_command_from_block(self, block, exclude_keys=None):
        """
        Build a subprocess-ready command list from a specific block.

        Parameters
        ----------
        block : dict
            One entry from self.blocks, representing a single executable command.
        exclude_keys : list or None
            Optional list of parameter keys to exclude.

        Returns
        -------
        list
            Command suitable for subprocess.run
        """

        if exclude_keys is None:
            exclude_keys = [
                'RADIAL', 'FalloffIsTrueSigma',
                'ActionIfGPUFails', 'UseGPU', 'FakeSIRTiterations'
            ]

        # program name comes from the block header
        program = block['header'].split()[0].lstrip('$')
        cmd = [program]

        for key in block['params']:
            if key in exclude_keys:
                continue

            val = block['params'][key]

            # --------------------------------------------------
            # Case 1: flag-only parameter (no argument)
            # --------------------------------------------------
            if val == '' or val is None or (isinstance(val, (list, tuple)) and len(val) == 0):
                cmd.append('-%s' % key)
                continue

            # --------------------------------------------------
            # Case 2: parameter with value
            # --------------------------------------------------
            sep = block['separators'].get(key)
            if sep is None:
                if isinstance(val, (list, tuple)):
                    sep = ','
                else:
                    sep = ''

            val_str = sep.join(self._val2str(val))

            cmd.append('-%s' % key)
            cmd.append(val_str)

        return cmd

    # ---------------------------
    # diagnostics / pretty print
    # ---------------------------
    def pretty_print_blocks(self):
        """
        Print a compact diagram of blocks and their keys.
        """
        out = []
        out.append('COMFILE: %s' % self.com_name)
        for i in range(len(self.blocks)):
            block = self.blocks[i]
            out.append('Block %d: %s' % (i, block['header']))
            keys = list(block['params'].keys())
            for k in keys:
                sep = block['separators'].get(k)
                out.append('  %s (sep=%s) -> %r' % (k, repr(sep), block['params'][k]))
            if block.get('footer'):
                out.append('  footer lines: %d' % len(block['footer']))
        print('\n'.join(out))
        return '\n'.join(out)
