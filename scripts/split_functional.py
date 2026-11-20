#!/usr/bin/env python3
"""One-shot generator: split src_orig/functional.f90 into the parent module
+ four F2008 implementation submodules (files already cut by region table
below; re-run only if the files are lost).

Design constraints validated on the densities split (commit 18793668):
  * procedures whose dummy/result declarations use explicit shapes bound to
    host variables (mv, nx, ...) stay in the parent: interface bodies cannot
    state such characteristics.
  * the parent's USE statements need only:-lists (gfortran 13 submodule
    attachment bug otherwise).
  * submodule files carry no top-level use lines; host association through
    the parent resolves everything.
  * Hephaestos substitutes the same keyword dictionary into every file
    (dispatch on name in src_heph/preprocess.py + multi-file substitute in
    heph_functional.ProcessFunctional).
"""
import re
import sys

SRC = 'src_orig/functional.f90'

# ---------------------------------------------------------------------- groups
GROUPS = {
    'functional_energy': [
        'CalcEnergy', 'CompSkyrme', 'PrintSkyrme',
        'CompKinetic_spwfs', 'CompKinetic_density',
        'calculate_epasta', 'calcspwfenergy',
        'CompCOMCorrection', 'calc_rotational_correction',
        'calc_vibrational_correction', 'inv_mat', 'classical_MOI',
        'CompStabilisingFactor', 'verify_COM_motion', 'calcElectronEnergy',
        'update_E_history',
    ],
    'functional_potentials': [
        'calcPotentials', 'calc_perturbed_potentials',
        'calc_perturbed_potentials_oneoff', 'combine_potentials',
        'add_constraint_potential', 'add_coulomb_potential',
        'add_cranking_potentials', 'precondition_potentials',
        'clean_potentials', 'save_potential_history',
    ],
    'functional_fields': [
        'Add_potentialvector', 'multiply_potentialvector', 'PVectorInproduct',
    ],
    'functional_io': [
        'WritePotentials', 'write_hdf5_potentials', 'ReadPotentials',
        'read_potentials_hdf5', 'hdf5_readpot', 'readpotentials_separate',
    ],
}
# Procedures that stay in the parent module (mv/nx-bound dummies or results)
PARENT_EXTRA = set()  # filled automatically below: all not in any group


def main():
    lines = open(SRC).read().splitlines()

    # ---- locate module structure anchors
    mod_start = next(i for i, l in enumerate(lines)
                     if re.match(r'module\s+functional', l, re.I))
    contains = next(i for i, l in enumerate(lines)
                    if re.match(r'^\s*contains\s*$', l, re.I))
    mod_end = next(i for i, l in enumerate(lines)
                   if re.match(r'end\s+module', l, re.I))
    spec = lines[mod_start:contains + 1]          # 'module' .. 'contains'

    # ---- locate every top-level procedure block
    starts = [i for i, l in enumerate(lines)
              if contains < i < mod_end
              and re.match(r'\s*(pure\s+)?(subroutine|function)\s+\w', l, re.I)]
    blocks = {}
    order = []
    for k, s in enumerate(starts):
        m = re.match(r'\s*(?:pure\s+)?(?:subroutine|function)\s+(\w+)', lines[s], re.I)
        name = m.group(1)
        decl = re.sub(r'\s+', ' ', lines[s]).strip()
        # capture multi-line statement up to the line ending without '&'
        j = s
        while lines[j].rstrip().endswith('&'):
            j += 1
        e = (starts[k + 1] if k + 1 < len(starts) else mod_end)
        blocks[name] = (s, j, e)
        order.append(name)

    grouped = set()
    for names in GROUPS.values():
        grouped.update(names)
    PARENT_EXTRA = [n for n in order if n not in grouped]

    # ---- generate interface bodies for grouped procedures
    def interface_of(name):
        s0, j0, e = blocks[name]
        stmt_join = []
        j = s0
        while True:
            stmt_join.append(lines[j].strip().rstrip('&').strip())
            if not lines[j].rstrip().endswith('&'):
                break
            j += 1
        stmt = ' '.join(stmt_join)
        dstr = ''
        m = re.search(r'\((.*)\)', stmt)
        if m:
            dstr = m.group(1).replace('&', '')
        # strip result(...) from dummy list handling
        rm = re.search(r'result\s*\(\s*(\w+)\s*\)', stmt, re.I)
        result_var = rm.group(1) if rm else None
        # dummy name extraction: cut at 'result(' if present
        dstr = re.sub(r'\bresult\s*\(.*', '', dstr)
        dummies = [d.strip() for d in dstr.split(',') if d.strip()]
        dummies = [d for d in dummies if '=' not in d] + \
                  [d.split('=')[0].strip() for d in dummies if d and '=' in d and not d.split('=')[1].strip().lower() in ('in', 'out', 'inout')]
        keep = [d for d in dummies if re.fullmatch(r'\w+', d)]
        if result_var:
            keep.append(result_var)
        # now scan the block spec for declaration lines mentioning dummies
        pde = re.compile(r'^(.*?)\b(' + '|'.join(keep) + r')\b')
        iface = []
        spec_end = None
        for i in range(j + 1, e):
            l = lines[i]
            if re.match(r'^\s*!', l) or not l.strip():
                continue
            if '::' in l and any(re.search(r'\b' + d + r'\b', l.split('::', 1)[1]) for d in keep):
                iface.append(l)
                # stop scanning decls once we leave spec part: heuristic =
                # first line after decls... we instead continue; executable
                # lines do not contain '::' with dummy names... risky but
                # compile-driven fixes follow.
        return stmt, iface, keep

    # ---- assertions on unmovable set
    unmovable = {'effmass_pot', 'set_coul', 'transfer_coulomb_mesh',
                 'INM_k2_pot_real', 'INM_k2_pot_complex',
                 'INM_k4_pot_real', 'INM_k4_pot_complex',
                 'apply_sphamil', 'delta_action', 'output_Edensity'}
    for u in unmovable:
        assert u not in grouped, f"{u} must stay in the parent (mv-bound)"

    # ---- write submodule files
    for fname, names in GROUPS.items():
        out = [f'submodule (functional) {fname}',
               '!==============================================================================',
               '!',
               f'! Implementation submodule of module functional: {fname}',
               '!  Interface block lives in src_orig/functional.f90.',
               '!',
               '!==============================================================================',
               'implicit none',
               'contains',
               '']
        for n in names:
            s, _, e = blocks[n]
            out.extend(lines[s:e])
            out.append('')
        out.append(f'end submodule {fname}')
        open(f'src_orig/{fname}.f90', 'w').write('\n'.join(out) + '\n')

    # ---- interface block for the parent
    iface_out = ['  interface',
                 '    !-------------------------------------------------------------------',
                 '    ! The following procedures are implemented in submodules (files',
                 '    ! src_orig/functional_{energy,potentials,fields,io}.f90).',
                 '    !-------------------------------------------------------------------']
    for gname, names in GROUPS.items():
        iface_out.append(f'    ! ---- {gname}')
        for n in names:
            stmt, decls, keep = interface_of(n)
            iface_out.append('    ' + stmt.replace(' &', '').strip())
            for d in decls:
                iface_out.append('    ' + d.strip())
            body = 'end function' if re.match(r'\s*(pure\s+)?function', stmt, re.I) else 'end subroutine'
            nm = re.match(r'\s*(?:pure\s+)?(?:subroutine|function)\s+(\w+)', stmt, re.I).group(1)
            iface_out.append(f'  {body} {nm}')
        iface_out.append('    ')
    iface_out.append('  end interface')
    iface_out.append('')

    # ---- parent file
    parent = (spec + iface_out)
    for n in PARENT_EXTRA:
        s, _, e = blocks[n]
        parent.extend(lines[s:e])
        parent.append('')
    parent.append('end module functional')
    open(SRC, 'w').write('\n'.join(parent) + '\n')

    print('wrote:', ', '.join(GROUPS), 'and parent', SRC)
    print('parent keeps:', ', '.join(PARENT_EXTRA))


if __name__ == '__main__':
    main()
