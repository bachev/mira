/*
 * Written by Thomas Pfisterer
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Thomas Pfisterer
 *
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 *
 * token definitions used in caf.files 
 * Version 0.91    27.07.99
 *         0.92    13.08.99  Umstellung auf Umgebungen
 *         0.93    15.09.99  einfache Erkennung von "non-ascii" files
 */





enum {
  token_error = 1,
  token_sequencename,
  token_identifier,
  token_quoted_text,
  token_number,
  token_type_read,
  token_type_contig,
  token_type_group,
  token_type_assembly,
  token_padded,
  token_unpadded,
  token_pstatus,  // 12
  token_state,
  token_asped,
  token_dye,
  token_dye_terminator,
  token_dye_primer,
  token_scf_file,
  token_primer,
  token_primer_unknown,
  token_primer_universal,
  token_primer_custom,  // 22
  token_template,
  token_insert_size,
  token_ligation,
  token_forward,
  token_reverse,
  token_seq_vector,
  token_sequencing_vector,
  token_clone_vector,
  token_clipping,
  token_golden_path,
  token_dna,
  token_sequence,
  token_quality,
  token_align_SCF,
  token_newline,
  token_assembled,
  token_staden_id,
  token_base_caller,
  token_tag,
  token_stolen,
  token_clone,
  token_ende
};

