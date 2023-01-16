# qpts
quantitatve analysis of gene ancestry

Q256005291344: tree_topology_v2.a.py first attempt of rewrtiting into functions
Q256007121136: v2.b - change algae dirs to dictionary
Q256007191654: v2.g - count_taxa() changed (not counting ',' but taxa)
Q256007232241: b = 101 addition
Q256010061643: aap_f.py - changes for Euglena gracilis Field
    Sorting of 18108 trees took 9110.300248597052 s.
Q256012061224: Script aap_f_first_hit.py renamed to aaq.py and tested on mariaDB data.
Q256101031442: Script aaq.py renamed to aar.py
    Trying to solve nearly-monophyly of EE taxa; the trees are counted multiple times, instead of once.
    Not started yet. 256101111656
Q256101171630: Script aar.py
    Annn replaced by Xnnn
Q256112191336 Script aar.py renamed to topology_tree_sorter_v4.0.py (changes for use in STEPA40db which usese 4codes!)
    improvements:
        1. deepcopy of ordereddict
        2. how to transform dict to 
Q256201031240 v4.1 establishing attempts - not working
Q256201041517 v4.2 establishing attempts
    print(categories, sep='\t') - vytiskne jako list, sep nemá žádnou funkci
    print(*categories, sep='\t') - vytiskne jako hodnoty oddělené '\t'
    troublesom tree:
        q1000034
        q1002049: 1 alga + 73 DErc#
Q256101141221 v4.21 not working for unknown reason, restored v4.2 from backup: 'Q:/backup/notepad++/tree_topology_sorter_v4.2.py.2019-01-09_214906.bak'
Q256101141246 v4.3 in case of sorting according to all phyla exception is raised (in tree q1000443) - proble caused by taxon bptx004183 (in count_target was wrongly the ignorecase option)
    working well
Q256201191436 v4.4 find_enter() returns the correct index
    tasks:
        1. The output of gene entrance must be in form neede by R, i.e.: [taxons] + [single None] + [internal nodes]
               the order of taxons is in the same order as in tree2lst(tree) output
               the order of nodes is in order from most outer to most inners
               between taxa and internal nodes must be one zero-element
Q256201231422 v4.5 - output gene entrance works well
    adding function count_ignored()
    it seems working well (tested on a single tree w/o DErc)
Q256201240706 v5.0 - adding plot option with RC
Q256203251359 v5.3 - seems working well - ran for all trees all categories and all euglenids, with bs=75; no tree, no proteins excluded (i.e. there could be some contamination and proteins without splice leader - esp. in Rapaza viridis - putative prey contamination)
    output: 20190325195344_13977_all_all_b=75
Q256203252001 'v5.31-RCo' - modification of 'unary' in R and of input tree (exclusion of Rhabdomonas costata)
    for formatting tip labels (italics, etc.) see
        https://guangchuangyu.github.io/software/ggtree/faq/
    It seems that the ignore_taxa problem is more complex; till now I've tried to skip non-ignored 'ignored taxa' results
        troublesome trees:
            q1006182
            q1009961
            q1014084
            q1022126
Q2562032271142 v5.32 checking if the results are exactly the same as in 5.3 (i.e. no ignored taxa)
    the result is exactly the same
Q256203281414 v5.33 concatenation of v5.32 and decontamination process of DERV and DEEE ('exclude_noSL_data_DERV_v0.2.py' and 'exclude_cont_data_DEEE_v0.3.py')
Q256203291409 v5.34
    1. check for DERV and DEEE decontamination; seems OK
    2. decontamination of MMETSP; works OK
Q256203311411 'v5.35' get rid of 102 and Qxxx.Qxxx (nonmonophyletic euglenids)
Q256204011559 'v5.36' correction for 102
Q256204171453 'v5.40' small changes to create output containing node of enter for each tree
Q256204171857 'v5.50' dealing with option for nested euglenids only; it seems I've tried already (around line 1008 in function 'analyse_sisterhood'
Q256204191855 'v5.51' added new result file contatining only nested
    one small issue, if inner branch has BS under treshold, it is excluded; not solving for JPD20019
        I think, it is not an issue: lower BS (-> multifurcation) in inner sisterhood could mean the target is not nested
    works probably well, but the results differ in also in target composition, which I did not expect
Q256204201219 'v5.52' creates also folder and copies only the nested tree into it; adds tree file size to the output, to easily choose small trees to check
    used for preselection of 'nested target trees'
Q256204201752 'v5.53' analyze nested groups only
    seems working good, surprisingly
    q1000054 - nested in nested with Kinetoplastea, cool
    q1003083 - target is not nested (in Kinetoplastea, but Kin is nested in target) :-(
        the same: q1009352 (Haptophyta)
    q1030817 - ochrophyta - wrong node? DEEE is nested in ochrophyta (denoted by group X001 - takes algal group), but the composition is for Q000 (takes node)
    q1000307 - two sister relationships are sister to each other (with Kinteoplastea) - nested, but strange
        the same: q1000058
Q256204202018 'v5.54' testing 4 wrongly sorted trees
    'remove_duplicite_line' should ignore bs, added to function
    q1003713 - wrong group
Q256204202317 'v5.55' - testing wrong trees
    Q000_absence must be global in all function
        after JPD2019 change to look into the target if Q000 presented
    preselection of 1427 trees
Q256204212046 'v5.56' - final for JPD2019
Q256205041758 'v5.57' - used for optimizing the BS treshold value
Q256205111747 'v6.00' - running for b = 50 and all 13977 trees (t = s)
Q256205122032 'v6.01' - adding mix - trees with mixed origin are copied into the folder mix
Q256205131758 'v6.02' - correction of one line ("results_nested_file = remove_duplicite_line(results_nested_file)" - in previous versions in else clausule, which is very wrong)
Q256205141724 'v6.03' - more clear view of results files, correct incorparation of mix in summary
Q256205301237 'v6.04' - first clear view for PSA 2019
    bs = 50
    nested relationship
    mix present, only algal groups sorted (x Metazoa as noise level?)
Q256205301908 'v6.05' count mix only if no specific algal group present for the Q_group in the tree
Q256205311302 'v6.06' header corrected
Q256205311307 'v6.07' 2072 input trees
Q256206011316 'v6.10' check of 6.05 novelty, seems to be wrong
q256206011734 'V6.11'

b = 111, i.e. the whole tree contain only taxa from one group of algae (except target)
not used: b = 101, i.e. tho whole tree contain only one taxon not belonging to the group of algae as the rest of the tree (except target).

#####################
### STEPA41db #######
#####################

Q256210161628 'v7.00' small changes 40db->41db (DEEZ -> DEEZ for Euglena longa)
Q256211041926 'v7.10' for unknown reason, there was used wrong file for summarization; corrected to unique row nested file in 7.01
    change order of algae
    change order of euglenids
Q256211070901 'v7.20' S/MeEGT is specified in the path to the script; if for nested/sister bootsraps
Q256301111412 'v7.30' run for 72263 out of 72854 for SUPP_Files_v2.0, b = 0
Q256301301450 'v7.30' run for 72263 out of 72854 for SUPP_Files_v2.0, b = 75 STA sister
Q256301301450 'v7.30' run for 72263 out of 72854 for SUPP_Files_v2.0, b = 50 STA nested
Q before 256308 - all subversion do not modify anything
Q256308111733 'v8.0' changed the sorting in nested analyses for mix category (previous version showed mix only if there was mix and mix (sister and sq(sister)), but in case, e.g. red, mix - it excluded the tree
    testing on six trees in case of 11-bins sorting and 5-bins (primary colours only) sorting
    tests:
                      |     S11SR75   | v7.30 | correct
            DEEE42318 | viridiplantae | nic   | mix
            DEEE42319 | OK            | OK    | viridiplantae 97, 56
            DEEE42681 | viridiplantae | nic   | mix
        should be solved by the condition: 'if prev_algal_group_is_mix or algal_group_is_mix:' 'label '2mix'
    Q256308170725 - it seems, the sorting is now OK
        variable copy_tree = True/False added to speed up the sorting:
            s11s75_1000_copy_tree_TRUE.out:  446.517438704 s
            s11s75_1000_copy_tree_FALSE.out: 431.124397157 s
                difference is 4 %
Q256312312225 'v8.1' small changes for STEPA42db
Q256509051521 'v8.11' renamed to 'edited_v8.11.py', adding prints to see hat is happening in the 2nd part of sorting, in which are some sorted results erroneously ommited (in the function 'remove_duplicite_line_mix()' [tag: RDLM]
Q256509051649 
    1. 'v8.2': 'tree_topology_sorter_v8.12_MTA_c.py' renamed to tts_v8.2.py,
    2. first try to correct partial blindfullness in RDLM
    3. first try to make general script which would be run by another scritp many times with all parameters
