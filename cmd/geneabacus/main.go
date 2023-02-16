//
// Copyright © 2015 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"

	"github.com/biogo/store/interval"

	"git.sr.ht/~vejnar/GeneAbacus/lib/esam"
	"git.sr.ht/~vejnar/GeneAbacus/lib/feature"
	"git.sr.ht/~vejnar/GeneAbacus/lib/profile"
)

var version = "DEV"

func parseStrand(strandRaw string) int8 {
	if strandRaw == "+" || strandRaw == "1" || strandRaw == "+1" {
		return 1
	}
	if strandRaw == "-" || strandRaw == "-1" {
		return -1
	}
	return 0
}

func main() {
	// Arguments: General
	var pathReport string
	var nWorker, verboseLevel int
	var appendOutput, verbose, printVersion bool
	flag.StringVar(&pathReport, "path_report", "", "Write report to path (stdout with -)")
	flag.IntVar(&nWorker, "num_worker", 1, "Number of worker(s)")
	flag.IntVar(&verboseLevel, "verbose_level", 0, "Verbose level")
	flag.BoolVar(&appendOutput, "append", false, "Append to output count and profile (default create)")
	flag.BoolVar(&verbose, "verbose", false, "Verbose")
	flag.BoolVar(&printVersion, "version", false, "Print version and quit")
	// Arguments: Input
	var pathSAMsRaw, pathBAMsRaw, rawSAMCmdIn, pathFeatures, formatFeatures, fonName, fonChrom, fonStrand, fonCoords, featureStrandRaw, pathFeaturesFilter, formatFeaturesFilter, fonNameFilter, fonChromFilter, fonStrandFilter, fonCoordsFilter, featureStrandRawFilter, libraryR1StrandRaw string
	var ignoreNHTag, paired, includeMissingInFilter bool
	flag.StringVar(&pathSAMsRaw, "path_sam", "", "Path to SAM file(s) (comma separated)")
	flag.StringVar(&pathBAMsRaw, "path_bam", "", "Path to BAM file(s) (comma separated)")
	flag.StringVar(&rawSAMCmdIn, "sam_command_in", "", "Command line to execute for opening each of the SAM file (comma separated)")
	flag.StringVar(&pathFeatures, "path_features", "", "Path to features file")
	flag.StringVar(&formatFeatures, "format_features", "FON", "Format of features file: 'FON' or 'tab'")
	flag.StringVar(&fonName, "fon_name", "transcript_stable_id", "FON key for feature name")
	flag.StringVar(&fonChrom, "fon_chrom", "chrom", "FON key for chromosome or locus")
	flag.StringVar(&fonStrand, "fon_strand", "strand", "FON key for strand")
	flag.StringVar(&fonCoords, "fon_coords", "exons", "FON key for coordinates (exons for example)")
	flag.StringVar(&featureStrandRaw, "feature_strand", "+", "Default feature strand (+ (+1) or - (-1))")
	flag.StringVar(&pathFeaturesFilter, "path_features_filter", "", "Path to features file (Filter)")
	flag.StringVar(&formatFeaturesFilter, "format_features_filter", "FON", "Format of features for Filter: 'FON' or 'tab'")
	flag.StringVar(&fonNameFilter, "fon_name_filter", "transcript_stable_id", "FON key for feature name for Filter")
	flag.StringVar(&fonChromFilter, "fon_chrom_filter", "chrom", "FON key for chromosome or locus for Filter")
	flag.StringVar(&fonStrandFilter, "fon_strand_filter", "strand", "FON key for strand for Filter")
	flag.StringVar(&fonCoordsFilter, "fon_coords_filter", "exons", "FON key for coordinates (exons for example) for Filter")
	flag.StringVar(&featureStrandRawFilter, "feature_strand_filter", "+", "Default feature strand for Filter (+ (+1) or - (-1))")
	flag.StringVar(&libraryR1StrandRaw, "read_strand", "", "Read 1 strand, i.e. + (+1) or - (-1) or unstranded if empty")
	flag.BoolVar(&paired, "paired", false, "Pair-end sequencing")
	flag.BoolVar(&ignoreNHTag, "ignore_nh_tag", false, "Ignore NH SAM tag and consider all alignment unique")
	flag.BoolVar(&includeMissingInFilter, "include_missing_in_filter", false, "Include missing feature in filter (present in main feature) as is")
	// Arguments: Read selection
	var minMappingQualityRaw, minOverlap, fragmentMinLength, fragmentMaxLength int
	var readLengthsRaw string
	var randProportionRaw float64
	var inProperPair bool
	flag.IntVar(&minMappingQualityRaw, "read_min_mapping_quality", 0, "Minimum read mapping quality")
	flag.IntVar(&minOverlap, "read_min_overlap", 10, "Minimum total overlap of the read with the feature interval(s)")
	flag.IntVar(&fragmentMinLength, "fragment_min_length", 0, "Minimum fragment length")
	flag.IntVar(&fragmentMaxLength, "fragment_max_length", 0, "Maximum fragment length")
	flag.StringVar(&readLengthsRaw, "read_length", "", "Read length(s) (comma separated)")
	flag.Float64Var(&randProportionRaw, "rand_proportion", -1., "Randomly select a proportion of all reads (from 0. to 1.)")
	flag.BoolVar(&inProperPair, "read_in_proper_pair", false, "Only read in proper pair (default: all pairs)")
	// Arguments: Counting
	var countPath, countMultisRaw, countTotalsRaw string
	var countTotalRealRead, countInProfile bool
	flag.StringVar(&countPath, "count_path", "counts.csv", "Path to counts output")
	flag.StringVar(&countMultisRaw, "count_multis", "1,2,900", "Read multiplicity to use for counting (comma separated)")
	flag.StringVar(&countTotalsRaw, "count_totals", "", "Totals (i.e. library size) for normalization (comma separated)")
	flag.BoolVar(&countTotalRealRead, "count_total_real_read", false, "Total read count is total number of read weighted (false) or not (true) by their multiplicity")
	flag.BoolVar(&countInProfile, "count_in_profile", false, "Only count reads included in the profile")
	// Arguments: Profiling
	var profilePathsRaw, profileTypeRaw, profileFormatsRaw string
	var profileMulti, profileOverhang, profileUntemplated, profileExtensionLength int
	var profilePositionFraction float64
	var profileNoUntemplated, profileNorm, profileNoCoordMapping bool
	flag.StringVar(&profilePathsRaw, "profile_paths", "profiles.bedgraph", "Path to profile output(s) (comma separated)")
	flag.StringVar(&profileTypeRaw, "profile_type", "", "Computing profiles of read distribution: 'first', 'last', 'first-last', 'position', 'all', 'all-extension' or 'all-slice'")
	flag.StringVar(&profileFormatsRaw, "profile_formats", "bedgraph", "Profile output format: 'bedgraph', 'binary' or 'csv' (comma separated)")
	flag.IntVar(&profileMulti, "profile_multi", 900, "Maximum alignment multiplicity to include a read in profile")
	flag.IntVar(&profileOverhang, "profile_overhang", 0, "Overhang length to add to each side of the profile")
	flag.IntVar(&profileUntemplated, "profile_untemplated", 0, "Remove max untemplated nucleotide")
	flag.IntVar(&profileExtensionLength, "profile_extension_length", 0, "Extension length for extension profile")
	flag.Float64Var(&profilePositionFraction, "profile_position_fraction", 0.5, "Fraction of position between start and end for position profile")
	flag.BoolVar(&profileNoUntemplated, "profile_no_untemplated", false, "Include only read w/o untemplated nucleotide in profile")
	flag.BoolVar(&profileNorm, "profile_norm", false, "Normalize profile counts with total reads")
	flag.BoolVar(&profileNoCoordMapping, "profile_no_coord_mapping", false, "Skip coordinate mapping from input to feature. Option specific to input and feature with the same coordinate system (e.g. genomic) only producing profile sense to the input. Used for genomic profile.")
	// Arguments: Output
	var pathMapping, pathSAMOutRaw string
	flag.StringVar(&pathMapping, "path_mapping", "", "Path to feature name(s) mapping (tabulated file)")
	flag.StringVar(&pathSAMOutRaw, "path_sam_out", "", "Path to output SAM file to save counted/mapped reads")
	// Arguments: Parse
	flag.Parse()

	// Version
	if printVersion {
		fmt.Println(version)
		os.Exit(0)
	}

	// Verbose
	if verbose && verboseLevel == 0 {
		verboseLevel = 1
	}

	// Max CPU
	runtime.GOMAXPROCS(nWorker * 2)

	// Time start
	var timeStart time.Time
	if verboseLevel > 0 {
		timeStart = time.Now()
	}

	// Check arguments
	if len(pathFeatures) == 0 {
		log.Fatal("No Feature input")
	} else if _, err := os.Stat(pathFeatures); os.IsNotExist(err) {
		log.Fatalln(pathFeatures, "not found")
	}

	// Parse raw arguments
	// pathSAMs
	var pathSAMs []esam.PathSAM
	var SAMCmdIn []string
	if len(pathSAMsRaw) > 0 {
		for _, p := range strings.Split(pathSAMsRaw, ",") {
			if _, err := os.Stat(p); os.IsNotExist(err) {
				log.Fatalln(p, "not found")
			} else {
				pathSAMs = append(pathSAMs, esam.PathSAM{Path: p, Binary: false})
			}
		}
		if len(rawSAMCmdIn) > 0 {
			SAMCmdIn = strings.Split(rawSAMCmdIn, ",")
		}
	}
	if len(pathBAMsRaw) > 0 {
		for _, p := range strings.Split(pathBAMsRaw, ",") {
			if _, err := os.Stat(p); os.IsNotExist(err) {
				log.Fatalln(p, "not found")
			} else {
				pathSAMs = append(pathSAMs, esam.PathSAM{Path: p, Binary: true})
			}
		}
	}
	if len(pathSAMs) == 0 {
		log.Fatal("No SAM/BAM input")
	}
	// libraryR1Strand
	libraryR1Strand := parseStrand(libraryR1StrandRaw)
	// readLengths
	var readLengths []int
	if len(readLengthsRaw) > 0 {
		for _, m := range strings.Split(readLengthsRaw, ",") {
			i, err := strconv.Atoi(m)
			if err != nil {
				log.Fatal(err)
			}
			readLengths = append(readLengths, i)
		}
		if verbose {
			timeNow := time.Now()
			fmt.Printf("%.1fmin - Keeping read of length: %v\n", timeNow.Sub(timeStart).Minutes(), readLengths)
		}
	}
	// fragmentLengths
	if fragmentMinLength > 0 && fragmentMaxLength > 0 {
		if verbose {
			timeNow := time.Now()
			fmt.Printf("%.1fmin - Keeping fragment from length %v to %v\n", timeNow.Sub(timeStart).Minutes(), fragmentMinLength, fragmentMaxLength)
		}
	} else if fragmentMinLength > 0 {
		if verbose {
			timeNow := time.Now()
			fmt.Printf("%.1fmin - Keeping fragment with min length %v\n", timeNow.Sub(timeStart).Minutes(), fragmentMinLength)
		}
	} else if fragmentMaxLength > 0 {
		if verbose {
			timeNow := time.Now()
			fmt.Printf("%.1fmin - Keeping fragment with max length %v\n", timeNow.Sub(timeStart).Minutes(), fragmentMaxLength)
		}
	}
	// minMappingQuality
	var minMappingQuality byte
	minMappingQuality = byte(minMappingQualityRaw)
	// randProportion
	var randProportion float32
	randProportion = float32(randProportionRaw)
	// countMultis
	var countMultis []int
	var profileMultiTotalCol int
	addProfileMulti := true
	for i, m := range strings.Split(countMultisRaw, ",") {
		im, err := strconv.Atoi(m)
		if err != nil {
			log.Fatal(err)
		}
		countMultis = append(countMultis, im)
		if im == profileMulti {
			addProfileMulti = false
			profileMultiTotalCol = i
		}
	}
	if addProfileMulti {
		countMultis = append(countMultis, profileMulti)
		profileMultiTotalCol = len(countMultis) - 1
	}
	profileMultiTotalCol = 1 + (2 * profileMultiTotalCol)
	// countTotals
	countTotals := make([]float64, 1+len(countMultis)*2)
	countTotalInput := false
	if countTotalsRaw != "" {
		for it, t := range strings.Split(countTotalsRaw, ",") {
			tf, err := strconv.ParseFloat(t, 64)
			if err != nil {
				log.Fatal(err)
			}
			countTotals[1+(2*it)] = tf
		}
		countTotalInput = true
	}
	// profileType
	var profileType int
	switch profileTypeRaw {
	case "first":
		profileType = profile.ProfileTypeFirst
	case "last":
		profileType = profile.ProfileTypeLast
	case "first-last":
		profileType = profile.ProfileTypeFirstLast
	case "position":
		profileType = profile.ProfileTypePosition
	case "all":
		profileType = profile.ProfileTypeAll
	case "all-slice":
		profileType = profile.ProfileTypeSplice
	case "all-extension":
		profileType = profile.ProfileTypeExtension
	default:
		profileType = profile.ProfileTypeNone
	}
	// Check arguments
	if (profileType == profile.ProfileTypeFirst || profileType == profile.ProfileTypeLast) && libraryR1Strand == 0 {
		log.Fatal("First and last position profile require stranded library (see read_strand option)")
	}
	// profilePaths
	var profilePaths []string
	profilePaths = strings.Split(profilePathsRaw, ",")
	// profileFormats
	var profileFormats []string
	profileFormats = strings.Split(profileFormatsRaw, ",")

	// Open features
	var features, featuresFilter, featuresMissing []feature.Feature
	var err error
	switch strings.ToLower(formatFeatures) {
	case "fon":
		features, err = feature.OpenFON(pathFeatures, fonName, fonChrom, fonStrand, fonCoords)
	case "tab":
		features, err = feature.OpenTAB(pathFeatures, parseStrand(featureStrandRaw))
	}
	if err != nil {
		log.Fatal(err)
	}

	// Open filter features
	var trees map[string]map[int8]*interval.IntTree
	if pathFeaturesFilter != "" {
		var featuresFilterRaw []feature.Feature
		var maxID uint32
		var err error
		switch strings.ToLower(formatFeaturesFilter) {
		case "fon":
			featuresFilterRaw, err = feature.OpenFON(pathFeaturesFilter, fonNameFilter, fonChromFilter, fonStrandFilter, fonCoordsFilter)
		case "tab":
			featuresFilterRaw, err = feature.OpenTAB(pathFeaturesFilter, parseStrand(featureStrandRawFilter))
		}
		if err != nil {
			log.Fatal(err)
		}
		// Check feature has corresponding filter-feature
		for _, feat := range features {
			found := false
			for _, featf := range featuresFilterRaw {
				if feat.Name == featf.Name {
					// Features used to filter and count/profile must have the same ID
					featf.ID = feat.ID
					// Keep max ID for features missing from filter
					if maxID < feat.ID {
						maxID = feat.ID
					}
					featuresFilter = append(featuresFilter, featf)
					found = true
					break
				}
			}
			if includeMissingInFilter && !found {
				fmt.Println("Warning: %s not found in filter, adding in filter as is", feat.Name)
				featuresMissing = append(featuresMissing, feat)
			}
		}
		// Add missing features only present in features but not in featuresFilter
		// Using maxID to put these features last
		for _, feat := range featuresMissing {
			maxID++
			feat.ID = maxID
			featuresFilter = append(featuresFilter, feat)
		}
		// Build feature trees
		trees, err = feature.BuildFeatTrees(featuresFilter)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		var err error
		// Build feature trees
		trees, err = feature.BuildFeatTrees(features)
		if err != nil {
			log.Fatal(err)
		}
	}

	// Open feature mapping
	var featuresMapping map[string]string
	if pathMapping != "" {
		featuresMapping, err = feature.OpenMapping(pathMapping)
	}
	if err != nil {
		log.Fatal(err)
	}

	// Output SAM
	var pathSAMOut esam.PathSAM
	if pathSAMOutRaw != "" {
		pathSAMOut = esam.PathSAM{Path: pathSAMOutRaw, Binary: false}
	}

	// Profile & Count alignments on Features
	nAlign, err := PConFeature(pathSAMs, SAMCmdIn, features, featuresMapping, trees, readLengths, fragmentMinLength, fragmentMaxLength, randProportion, paired, libraryR1Strand, ignoreNHTag, inProperPair, minMappingQuality, minOverlap, countMultis, countTotals, countTotalInput, countTotalRealRead, countInProfile, countPath, profileType, profileMulti, profileOverhang, profileNoCoordMapping, profileUntemplated, profileNoUntemplated, profileExtensionLength, profilePositionFraction, profileNorm, profileMultiTotalCol, profilePaths, profileFormats, appendOutput, pathReport, pathSAMOut, nWorker, timeStart, verboseLevel)
	if err != nil {
		log.Fatal(err)
	}

	// Verbose
	if verboseLevel > 0 {
		timeEnd := time.Now()
		fmt.Printf("%.1fmin - Done %d align.\n", timeEnd.Sub(timeStart).Minutes(), nAlign)
	}
}
