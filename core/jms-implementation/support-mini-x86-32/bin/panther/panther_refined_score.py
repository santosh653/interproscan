#parse the hmmer3 tblout file
#run the pfsearch binary
#gift nuka
#jan 2015, mar 2015


import subprocess
import sys
import os.path
import os
import re
import time

import getopt

from tempfile import NamedTemporaryFile

#total runtime for hmmer runs
total_run_time = 0

#reformat pfsearch alignment output to appear on one line
def clean_output(lines):
  line_list = iter(lines.splitlines())
  formatted_lines  = ''

  for line in line_list:
    if not line.strip():
      continue
    if 'match_nb' in line:
      line = re.sub('match_nb=\S+\s+match_type=\S+\s+', '', line)
      if 'match_nb' in line:
        #print 'regex failed ... ', line
        line = re.sub('match_nb=\S+\s+', '', line)
        line = re.sub('match_type=\S+\s+', '', line)
        #print 'regex may hav failed, print line again: ', line
	#else print 'regex worked ... ', line
    if formatted_lines and len(line) > 0 and line[0] == '>':
      formatted_lines += '\n'
      formatted_lines += line
    else:
      formatted_lines += line
  #print formatted_lines
  if not formatted_lines.strip():
    formatted_lines = ''
  else:
    formatted_lines += '\n'
  #print formatted_lines
  return formatted_lines

##get the hamap profiles from the
def get_hamap_profile(profiles_list_filename):
    profiles = {}
    lines = []
    if not os.path.isfile(profiles_list_filename):
        return profiles

    temp_err_file = profiles_list_filename + '-filtered'
    with open(profiles_list_filename, "r") as profile_list:
        for line in profile_list:
            line = line.strip()
            if not line.startswith('#'):
                lines.append(line)
                m = re.search('^(\S+)\s+(\S+)\s+(\S+)\s+(.*)', line)
                hit_line = ''
                if m:
                    seq_id = m.group(1)
                    profile = m.group(3)
                    profile_path = model_dir + '/' + profile + ".prf"
                    #print line
                    hit_line = 'ok - ' + line
                    if profile in profiles:
                       profiles[profile].extend([seq_id])
                    else:
                       profiles[profile] = [profile_path, seq_id]
                       #print profiles[profile]
                else:
                    print('something wrong ' + line)
                    sys.stderr.write('something wrong ' + line)
                    hit_line = 'not ok - ' + line
                #append_to_file(temp_err_file, hit_line+'\n')
    key_count = len(list(profiles.keys()))
    append_to_file(temp_err_file, 'profile hits: ' + str(key_count) + '\n')
    return profiles

def get_sequences_from_fasta(fasta_file):
    line_count = 0
    print ('get_sequences_from_fasta')
    seq_id ='000'
    with open(fasta_file, 'r') as fasta:
        fasta_dict = {}
        for line in fasta:
            line = line.strip()
            if line == '':
                continue
            if line.startswith('>'):
                seq_id = line.lstrip('>')
                seq_id = re.sub('\..*', '', seq_id)
                fasta_dict[seq_id] = ''
            else:
                #print('seq_id: ' + seq_id)
                fasta_dict[seq_id] += line + '\n'
                line_count = line_count + 1
                #if len(line) > 80:
                #    raise ValueError('Input fasta file format problem for hmmer, line length greater than 80 ')
    print('done get_sequences_from_fasta')
    return fasta_dict

def get_sequences_for_profile(key_list, seqs_dict):
    sequences = ''
    #print len(key_list),':', key_list
    for key in key_list:
        if key in seqs_dict:
            #print "key found : ", key
            value = seqs_dict[key]
            sequences += '>' + key + '\n' + value
            if not value.endswith('\n'):
                sequences += '\n'
	    #print 'ok', sequences
        else:
            print("key not found : " + key)

    #print 'ok', sequences
    return sequences

def get_sequence(key, seqs_dict):
    sequences = ''
    if key in seqs_dict:
       value = seqs_dict[key]
       sequences += '>' + key + '\n' + value
    else:
       print("key not found : " + key)
    return sequences

def create_temp_file(filename, temp_dir):
    file_prefix = filename + '.'
    #temp_dir = '/nfs/nobackup/interpro/nuka/i5/build/release-5.18-hamap/may06-test/temp/tmp'
    f = NamedTemporaryFile(delete=False, prefix=file_prefix, dir=temp_dir)
    return f.name

def append_to_file(filename, output):
    with open(filename, 'a') as out_file:
        out_file.write(output)

def write_to_file(filename, output):
    with open(filename, 'w') as seq_file:
        seq_file.write(output)


def get_query_name(hmma):
    """
      get the panther family name from the query target
    """
    hmma_list = hmma.split ('.')
    if len(hmma_list) > 2:
        hmm_fam  = hmma_list[0]
        hmm_sf = hmma_list[1]
        something_else = hmma_list[2]
    elif len(hmma_list) == 2:
        hmm_fam  = hmma_list[0]
        something_else = hmma_list[1]
    hmm_id = hmm_fam
    if hmm_sf and hmm_sf.startswith("SF"):
        hmm_id = hmm_fam + ':' + hmm_sf
    return hmm_id

def get_panther_families(names_tab_file):
    panther_families = {}
    with open(names_tab_file, 'r') as infile:
        for line in infile:
            hmm_fam, family_name = line.split('\t')
            hmm_id = get_query_name(hmm_fam)
            panther_families[hmm_id] = family_name.strip()
    return panther_families

def get_singleton_panther_families(family_names_file):
    singleton_panther_families = []
    with open(family_names_file, 'r') as infile:
        for line in infile:
            family_name = line.strip()
            if family_name:
                singleton_panther_families.append(family_name)
    return singleton_panther_families

def append_to_match_list(all_scores, seq_id, item):
    updated_raw_matches = []
    if seq_id in all_scores:
        raw_matches = all_scores[seq_id]
        if  not raw_matches:
            #wwe expect hits for each seqid
            print_error ('some problem with the hmmer output, ...')
            sys.exit(4)
        raw_matches.append(item)
        updated_raw_matches = raw_matches
    else:
        updated_raw_matches = [item]
    return updated_raw_matches

#return true if two int ranges overlap
def location_overlaps(location1, location2):
    start_1,end_1 = list(map(int, location1.split("-")))
    start_2,end_2 = list(map(int,location2.split("-")))
    return max(start_1,start_2) <= min(end_1,end_2)

def get_match_groups(best_hits):
    group_hits = {}
    for el in best_hits:
        hmm_id = el[0]
        if hmm_id in group_hits:
            el_hits = group_hits[hmm_id]
            if el not in el_hits:
                el_hits.append(el)
        else:
            group_hits[hmm_id] = [el]
    return group_hits

def get_family_matches(matches):
    family_matches = []
    for match in matches:
        hmm_id = match[0]
        hmm_fam = hmm_id
        evalue = match[2]
        if 'SF' in  hmm_id:
            hmm_fam, hmm_sf = hmm_id.split(':')
        new_match = [hmm_fam,evalue]
        family_matches.append(new_match)
    return family_matches


def get_filtered_best_hits(best_hits):
    # hmm_hit = [hmm_id, description, float(eVal), float(score), location]
    filtered_best_hits = []
    all_best_hits = best_hits
    hits_per_hmmid = get_match_groups(best_hits)
    if len(list(hits_per_hmmid.keys())) == 1:
        return best_hits
    for hmmid in hits_per_hmmid:
        hmmid_hits =  hits_per_hmmid[hmmid]
        if ":SF" in hmmid:
            filtered_best_hits.extend(hmmid_hits)
            continue
        for el in hmmid_hits:
            base_el_location = el[4]
            overlaps = False
            for other_el in all_best_hits:
                if el == other_el:
                    continue
                if el[0] == other_el[0]:
                    continue
                if location_overlaps(base_el_location, other_el[4]):
                    overlaps = True
                    break
            if overlaps:
                all_best_hits.remove(el)
            else:
                filtered_best_hits.append(el)
    return filtered_best_hits

def get_best_hits(matches, evalue_cutoff):
    """
      get the best hit: one with smallest evalue and highest score.
      if at least two hits have same evalue and score, report them all
    """
    best_hits = []
    evalue_sorted = sorted(matches, key=lambda x: x[2])
    print(len(evalue_sorted))
#     for el in evalue_sorted:
#       print (el)
#       break
    best_evalue = evalue_sorted[0][2]
    if best_evalue <= evalue_cutoff:
        scores = []
        for el in evalue_sorted:
            if len(el) >= 4:
                if best_evalue == el[2]:
                    scores.append(el)
        score_sorted = sorted(scores, key=lambda x: x[3], reverse=True)
        best_score = score_sorted[0][3]
        for el in score_sorted:
            if len(el) >= 4:
                if best_score == el[3]:
                    best_hits.append(el)
    return best_hits

def parse_domtblout(domtblout, panther_families, run_mode):
    all_scores = {}
    with open(domtblout, 'r') as infile:
        for line in infile:
            if not line.strip():
                #possible empty lines
                continue
            if not line.startswith('#'):
                if run_mode == 'hmmscan':
                    hmma, acc2, qlen, seqid, acc1, tlen, eVal, score, bias, num, num_t, cEval, iEval, dscore, dbias, hmm_f, hmm_t, ali_f, ali_t, env_f, env_t, accuracy, desc = line.split()
                elif run_mode == 'hmmsearch':
                    seqid, acc1, tlen, hmma, acc2, qlen, eVal, score, bias, num, num_t, cEval, iEval, dscore, dbias, hmm_f, hmm_t, ali_f, ali_t, env_f, env_t, accuracy, desc = line.split()
                else:
                    print_error ("error: run_mode is invalid: " + run_mode)
                hmm_id = get_query_name(hmma)
                description = panther_families[hmm_id]
                hmm_location = str(hmm_f) + '-' + str(hmm_t)
                ali_location = str(ali_f) + '-' + str(ali_t)
                env_location = str(env_f) + '-' + str(env_t)
                hmm_hit = [hmm_id, description, float(eVal), float(score), hmm_location, ali_location, env_location, qlen]
                matches = append_to_match_list(all_scores, seqid, hmm_hit)
                all_scores[seqid] = matches
    return all_scores


def run_hmmer_binary(hmmer_cmd, temp_dir, hmmer_options, panther_family_hmm_path, sequence_key, input_fasta_sequences, panther_families):
    global total_run_time
    space = ' '
    temp_file = temp_dir + '/' + sequence_key
    domtbl_file = temp_file + '.domtbl.out'
    temp_output_file = temp_file + '.raw.out'
    write_to_file(temp_file, input_fasta_sequences)
    panther_command = hmmer_cmd + space + hmmer_options + ' -o ' + temp_output_file + ' --domtblout ' + domtbl_file +  space + panther_family_hmm_path + space + temp_file
    print (panther_command)
    start_time = time.time()
    outflag = os.system(panther_command)
    end_time = time.time()
    run_time = end_time - start_time
    print ('run hmmscan for :' + str(run_time))
    total_run_time = total_run_time + run_time
    run_mode = 'hmmscan'
    mini_scores = parse_domtblout(domtbl_file, panther_families, run_mode)
    return mini_scores

def old_pfsearch_binary_run():
    count = 0
    temp_file_list = []
    stats_filename = input_fasta_file + ".stats"
    #base_job_dir = os.path.dirname(stats_filename)
    temp_dir = input_fasta_file + '-tmp'
    os.makedirs(temp_dir)
    #prf_half = int(len(profiles) / 2)
    temp_err_file = profiles_list_filename + '-filtered'

    append_to_file(temp_err_file, 'profile hits in run_pfsearch: ' + str(len(list(profiles.keys()))) + '\n')

    append_to_file(temp_err_file, 'sequence count in run_pfsearch: ' + str(len(list(seqs_dict.keys()))) + '\n')
    get_seq_time = 0
    keys = list(profiles.keys())
    write_to_file(temp_err_file + '-keys', 'profile keys in run_pfsearch: ' + str(len(keys)) + '\n')
    for prf in profiles:
        prf_seqs = ' '.join(profiles[prf])
        prf_seqs_count = len(profiles[prf][1:])
        prf_out = 'processing profile #:' + str(count) + ' ' +  prf + ' seqs:' + str(prf_seqs_count) + ' - ' + prf_seqs + '  \n'
        append_to_file(temp_err_file, prf_out)
        if prf == 'MF_00005':
            #append_to_file(temp_err_file, 'processing MF_00005 \n')
            mf_output = ' '.join(profiles[prf])
            #append_to_file(temp_err_file, mf_output)
            #append_to_file(temp_err_file, ' ok \n ')
        sequence_ids = profiles[prf][1:]
        get_seq_start_time = time.time()
        input_fasta_sequences = get_sequences_for_profile(sequence_ids, seqs_dict)
        #print input_fasta_sequences
        #temp_file = create_temp_file(prf, temp_dir)
        temp_file = temp_dir + '/' + prf
        temp_output_file = temp_file + '.raw.out'
        #print 'temp file is ', temp_file
        write_to_file(temp_file, input_fasta_sequences)
        get_seq_end_time = time.time()
        time_to_get_seqences = get_seq_end_time - get_seq_start_time
        #print 'time to get ', len(sequence_ids), 'sequences for ', prf, ': ', time_to_get_seqences
        get_seq_time += time_to_get_seqences
        temp_file_list.append(temp_file)

        comd_to_run = []
        comd_to_run = arg_list[command_index:]
        #comd_to_run.extend([profiles[prf][0],fasta_file])
        comd_to_run.extend([profiles[prf][0],temp_file])
        cmd_string = ' '.join(comd_to_run)
        #print "command to run: ",  cmd_string
        count = count  + 1
        #append_to_file(temp_err_file, "command to run: " +  cmd_string + ' \n')
        if not os.path.isfile(profiles[prf][0]):
            #profile not available
            continue
        output = subprocess.check_output(comd_to_run, universal_newlines=True)
        #append_to_file(temp_err_file, "command to run: " +  cmd_string + ' \n')
        #if prf == 'MF_00005':
        #    append_to_file(temp_err_file, "command to run: " +  cmd_string + ' \n')
        #    with open(temp_output_file, 'a') as out_temp_file:
        #        out_temp_file.write(output)
        output = clean_output(output)
        if output.strip():
            with open(output_file, 'a') as out_file:
                out_file.write(output)
        #if count > prf_half:
        #    break

    for tempfile in temp_file_list:
        #print tempfile
        #todo remove comment
        #os.unlink(tempfile)
        testcount = 1

    append_to_file(temp_err_file, 'completed running thru ' + str(count) + ' profiles \n')
    with open(stats_filename, 'w') as stats_file:
        stats_file.write('Total time to get and write ' + str(len(temp_file_list)) + ' seq files :' + str(get_seq_time * 1000) + " ms \n")
    return count

def print_list(seq_key, best_hits):
    out_str_list = []
    for el in best_hits:
        out_str_list.append (seq_key + "\t" + "\t".join(str(x) for x in el) )
        #promote parent of subfamily
        if ':SF' in el[0]:
            el[0] = el[0].split (':')[0]
            out_str_list.append ( seq_key + "\t" + "\t".join(str(x) for x in el) )
    return out_str_list

def usage():
    sys.stderr.write("usage: python panther_score.py -d domtbl.out -m hmmscan|hmmsearch -n names.tab -e 0.001\n" )


if __name__ == "__main__":

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:d:i:n:e:m:o:p:s:c:", ["help", "input-fasta","domtbl", "hmmscan-cmd", "runmode", "evalue-cutoff","output", "panther-models-dir", "singleton-fam"])
    except getopt.GetoptError as err:
        print(err)  # should print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for option, arg_value in opts:
        if option in ("-d", "--domtbl"):
            domtblout = arg_value
        elif option in ("-h", "--help"):
            usage()
            sys.exit()
        elif option in ("-i", "--input-fasta"):
            fasta_file = arg_value
        elif option in ("-m", "--runmode"):
            run_mode = arg_value
        elif option in ("-c", "--hmmscan-cmd"):
                    hmmscan_cmd = arg_value
        elif option in ("-p", "--panther-models-dir"):
            panther_models_dir = arg_value
        elif option in ("-n", "--nametab"):
            names_tab_file = arg_value
        elif option in ("-s", "--singleton-fam"):
            singleton_familes = arg_value
        elif option in ("-e", "--evalue-cutoff"):
            evalue_cutoff = float(arg_value)
        elif option in ("-o", "--output"):
            output_file = arg_value
        else:
            assert False, "unhandled option"

    try:
        if not (domtblout and run_mode and names_tab_file and output_file):
            print_error("provide expected options")
            usage()
            sys.exit(3)
        if not evalue_cutoff:
            evalue_cutoff = float(1e-11)
    except NameError:
        print ("provide the required parameters")
        usage()
        sys.exit(3)

    try:
        start_time = time.time()

        cmd_str = ' '.join(sys.argv)
        print (cmd_str)
        panther_families = get_panther_families(names_tab_file)
        singleton_models =  get_singleton_panther_families(singleton_familes)
        all_scores = parse_domtblout(domtblout, panther_families, run_mode)

        #create the output file in case we don't have any matches
        open(output_file, 'w').close()
        #get the protein sequences
        print("now process the refined details")
        seqs_dict = get_sequences_from_fasta(fasta_file)
        print(seqs_dict.keys())

        #panther_models_hit = get_panther_models(initial_hmmer_output)
        end_time = time.time()
        read_file_time = end_time - start_time
        start_time = time.time()
        temp_dir = 'temp'
        #hmmer_cmd = '/home/gift/ebi/i5/test/may21/bin/hmmer/hmmer3/3.1b1/hmmscan'
        #hmmscan_cmd
        hmmer_run_count = 0
        sngl_fam_matches = 0
        hmmer_options = '-E 0.001 --domE 0.00000001 --incdomE 0.00000001 -Z 65000000 --cpu 1'
        with open(output_file, 'a') as outf:
            for seq_key in all_scores.keys():
                #stats_filename = fasta_file + ".stats"
                #run the pfsearch binary
                matches = all_scores[seq_key]
                #best_hits = get_best_hits(matches,  evalue_cutoff)
                #for sub families only include families
                family_matches = get_family_matches(matches)
                ## maybe remove duplicates
                processed_hmm = []
                seq_scores = []
                for match in matches:
                    hmm_id = match[0]
                    if hmm_id in processed_hmm:
                        continue
                    panther_family_hmm_path = panther_models_dir + '/' + hmm_id + '.sf.hmm'
                    #evalue = match[2]
                    #if evalue <= evalue_cutoff: TODO
                    if hmm_id in singleton_models:
                        seq_scores.append(match)
                        sngl_fam_matches += 1
                        continue
                    else:
                        processed_hmm.append(hmm_id)
                        #process this id
                        #temp_seq_key = 'UPI0002E0D40B'
                        input_fasta_sequences = get_sequence(seq_key, seqs_dict)
                        mini_scores = run_hmmer_binary(hmmscan_cmd, temp_dir, hmmer_options, panther_family_hmm_path, seq_key, input_fasta_sequences, panther_families)
                        seq_scores.extend(mini_scores[seq_key])
                        hmmer_run_count += 1
                if len(seq_scores) > 0:
                    seq_best_hits = get_best_hits(seq_scores,  evalue_cutoff)
                    seq_filtered_best_hits = get_filtered_best_hits(seq_best_hits)
                    seq_out_str_list = print_list(seq_key, seq_filtered_best_hits)
                    for out_str in seq_out_str_list:
                        outf.write(out_str + '\n')
                    #print('finished processing .. sequence: ' + str(seq_key))
        #sys.stdout.flush()
        print('Total hmmerscan runs: ' + str(hmmer_run_count))
        print('Total hmmerscan runtime: ' + str(total_run_time))
        print('Avge hmmerscan runtime: ' + str(total_run_time/hmmer_run_count))
        print('sngl_fam_matches: ' + str(sngl_fam_matches))
    except:
        print(sys.version)
        print("Unexpected error: ")
        print(sys.exc_info())