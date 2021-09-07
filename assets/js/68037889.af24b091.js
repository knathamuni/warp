"use strict";(self.webpackChunkwebsite_2=self.webpackChunkwebsite_2||[]).push([[1798],{3905:function(e,t,a){a.d(t,{Zo:function(){return m},kt:function(){return c}});var n=a(7294);function i(e,t,a){return t in e?Object.defineProperty(e,t,{value:a,enumerable:!0,configurable:!0,writable:!0}):e[t]=a,e}function l(e,t){var a=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),a.push.apply(a,n)}return a}function r(e){for(var t=1;t<arguments.length;t++){var a=null!=arguments[t]?arguments[t]:{};t%2?l(Object(a),!0).forEach((function(t){i(e,t,a[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(a)):l(Object(a)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(a,t))}))}return e}function o(e,t){if(null==e)return{};var a,n,i=function(e,t){if(null==e)return{};var a,n,i={},l=Object.keys(e);for(n=0;n<l.length;n++)a=l[n],t.indexOf(a)>=0||(i[a]=e[a]);return i}(e,t);if(Object.getOwnPropertySymbols){var l=Object.getOwnPropertySymbols(e);for(n=0;n<l.length;n++)a=l[n],t.indexOf(a)>=0||Object.prototype.propertyIsEnumerable.call(e,a)&&(i[a]=e[a])}return i}var s=n.createContext({}),p=function(e){var t=n.useContext(s),a=t;return e&&(a="function"==typeof e?e(t):r(r({},t),e)),a},m=function(e){var t=p(e.components);return n.createElement(s.Provider,{value:t},e.children)},u={inlineCode:"code",wrapper:function(e){var t=e.children;return n.createElement(n.Fragment,{},t)}},d=n.forwardRef((function(e,t){var a=e.components,i=e.mdxType,l=e.originalType,s=e.parentName,m=o(e,["components","mdxType","originalType","parentName"]),d=p(a),c=i,f=d["".concat(s,".").concat(c)]||d[c]||u[c]||l;return a?n.createElement(f,r(r({ref:t},m),{},{components:a})):n.createElement(f,r({ref:t},m))}));function c(e,t){var a=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var l=a.length,r=new Array(l);r[0]=d;var o={};for(var s in t)hasOwnProperty.call(t,s)&&(o[s]=t[s]);o.originalType=e,o.mdxType="string"==typeof e?e:i,r[1]=o;for(var p=2;p<l;p++)r[p]=a[p];return n.createElement.apply(null,r)}return n.createElement.apply(null,a)}d.displayName="MDXCreateElement"},2305:function(e,t,a){a.r(t),a.d(t,{frontMatter:function(){return r},contentTitle:function(){return o},metadata:function(){return s},toc:function(){return p},default:function(){return u}});var n=a(7462),i=a(3366),l=(a(7294),a(3905)),r={sidebar_position:2},o="VCF Overview: Illumina Genotyping Array",s={unversionedId:"Pipelines/Illumina_Genotyping_Arrays_Pipeline/Illumina_genotyping_array_spec",id:"Pipelines/Illumina_Genotyping_Arrays_Pipeline/Illumina_genotyping_array_spec",isDocsHomePage:!1,title:"VCF Overview: Illumina Genotyping Array",description:"The Illumina Genotyping Array Pipeline v1.11.0 pipeline produces a VCF (Variant Call Format) output with data processing and sample-specific genotype information. The VCF follows the format listed in the VCF 4.2 specification, but additionally contains fields and attributes that are unique to the Arrays pipeline.",source:"@site/docs/Pipelines/Illumina_Genotyping_Arrays_Pipeline/Illumina_genotyping_array_spec.md",sourceDirName:"Pipelines/Illumina_Genotyping_Arrays_Pipeline",slug:"/Pipelines/Illumina_Genotyping_Arrays_Pipeline/Illumina_genotyping_array_spec",permalink:"/warp/docs/Pipelines/Illumina_Genotyping_Arrays_Pipeline/Illumina_genotyping_array_spec",editUrl:"https://github.com/broadinstitute/warp/edit/develop/website/docs/Pipelines/Illumina_Genotyping_Arrays_Pipeline/Illumina_genotyping_array_spec.md",tags:[],version:"current",lastUpdatedBy:"Farzaneh Khajouei",lastUpdatedAt:1631050942,formattedLastUpdatedAt:"9/7/2021",sidebarPosition:2,frontMatter:{sidebar_position:2},sidebar:"docsSidebar",previous:{title:"Illumina Genotyping Array Overview",permalink:"/warp/docs/Pipelines/Illumina_Genotyping_Arrays_Pipeline/IlluminaGenotypingArray.documentation"},next:{title:"Optimus Overview",permalink:"/warp/docs/Pipelines/Optimus_Pipeline/README"}},p=[{value:"Meta information fields",id:"meta-information-fields",children:[]},{value:"Header line columns",id:"header-line-columns",children:[{value:"FILTER",id:"filter",children:[]},{value:"FORMAT (genotype)",id:"format-genotype",children:[]},{value:"INFO",id:"info",children:[]}]}],m={toc:p};function u(e){var t=e.components,a=(0,i.Z)(e,["components"]);return(0,l.kt)("wrapper",(0,n.Z)({},m,a,{components:t,mdxType:"MDXLayout"}),(0,l.kt)("h1",{id:"vcf-overview-illumina-genotyping-array"},"VCF Overview: Illumina Genotyping Array"),(0,l.kt)("p",null,"The ",(0,l.kt)("a",{parentName:"p",href:"https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/genotyping/illumina/IlluminaGenotypingArray.wdl"},"Illumina Genotyping Array Pipeline")," v1.11.0 pipeline produces a VCF (Variant Call Format) output with data processing and sample-specific genotype information. The VCF follows the format listed in the ",(0,l.kt)("a",{parentName:"p",href:"https://samtools.github.io/hts-specs/VCFv4.2.pdf"},"VCF 4.2 specification"),", but additionally contains fields and attributes that are unique to the Arrays pipeline."),(0,l.kt)("p",null,"This document describes the Array pipeline\u2019s unique VCF fields and attributes that are not listed in the standard VCF specification. To learn more about the pipeline, see the ",(0,l.kt)("a",{parentName:"p",href:"/warp/docs/Pipelines/Illumina_Genotyping_Arrays_Pipeline/IlluminaGenotypingArray.documentation"},"Illumina Genotyping Array Pipeline Overview"),"."),(0,l.kt)("div",{className:"admonition admonition-tip alert alert--success"},(0,l.kt)("div",{parentName:"div",className:"admonition-heading"},(0,l.kt)("h5",{parentName:"div"},(0,l.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,l.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"12",height:"16",viewBox:"0 0 12 16"},(0,l.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.5 0C3.48 0 1 2.19 1 5c0 .92.55 2.25 1 3 1.34 2.25 1.78 2.78 2 4v1h5v-1c.22-1.22.66-1.75 2-4 .45-.75 1-2.08 1-3 0-2.81-2.48-5-5.5-5zm3.64 7.48c-.25.44-.47.8-.67 1.11-.86 1.41-1.25 2.06-1.45 3.23-.02.05-.02.11-.02.17H5c0-.06 0-.13-.02-.17-.2-1.17-.59-1.83-1.45-3.23-.2-.31-.42-.67-.67-1.11C2.44 6.78 2 5.65 2 5c0-2.2 2.02-4 4.5-4 1.22 0 2.36.42 3.22 1.19C10.55 2.94 11 3.94 11 5c0 .66-.44 1.78-.86 2.48zM4 14h5c-.23 1.14-1.3 2-2.5 2s-2.27-.86-2.5-2z"}))),"How do I view a VCF file?")),(0,l.kt)("div",{parentName:"div",className:"admonition-content"},(0,l.kt)("p",{parentName:"div"},"A VCF is a tab delimited text file that you can open with either a text editor or with a spreadsheet editor (like Microsoft Excel)."))),(0,l.kt)("h2",{id:"meta-information-fields"},"Meta information fields"),(0,l.kt)("p",null,"Each VCF has meta information fields with attributes that generally describe the processing of the sample in the array."),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"arrayType - Name of the genotyping array that was processed."),(0,l.kt)("li",{parentName:"ul"},"autocallDate - Date that the genotyping array was processed by \u2018autocall\u2019 (aka gencall), the Illumina genotype calling software"),(0,l.kt)("li",{parentName:"ul"},"autocallGender - Gender (sex) that autocall determined for the sample processed"),(0,l.kt)("li",{parentName:"ul"},"autocallVersion - Version of the autocall/gencall software used"),(0,l.kt)("li",{parentName:"ul"},"chipWellBarcode - Chip well barcode (a unique identifier for sample as processed on a specific location on the Illumina genotyping array)"),(0,l.kt)("li",{parentName:"ul"},"clusterFile - Cluster file used"),(0,l.kt)("li",{parentName:"ul"},"extendedIlluminaManifestVersion - Version of the \u2018extended Illumina manifest\u2019 used by the VCF - generation software."),(0,l.kt)("li",{parentName:"ul"},"extendedManifestFile - File name of the \u2018extended Illumina manifest\u2019 used by the VCF generation software"),(0,l.kt)("li",{parentName:"ul"},"fingerprintGender - Gender (sex) determined using an orthogonal fingerprinting technology, populated by an optional parameter used by the VCF generation software"),(0,l.kt)("li",{parentName:"ul"},"gtcCallRate - GTC call rate of the sample processed that is generated by the autocall/gencall software and represents the fraction of callable loci that had valid calls"),(0,l.kt)("li",{parentName:"ul"},"imagingDate - Creation date for the chip well barcode IDATs (raw image scans)"),(0,l.kt)("li",{parentName:"ul"},"manifestFile - Name of the Illumina manifest (.bpm) file used by the VCF generation software"),(0,l.kt)("li",{parentName:"ul"},"sampleAlias - Sample name")),(0,l.kt)("div",{className:"admonition admonition-tip alert alert--success"},(0,l.kt)("div",{parentName:"div",className:"admonition-heading"},(0,l.kt)("h5",{parentName:"div"},(0,l.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,l.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"12",height:"16",viewBox:"0 0 12 16"},(0,l.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.5 0C3.48 0 1 2.19 1 5c0 .92.55 2.25 1 3 1.34 2.25 1.78 2.78 2 4v1h5v-1c.22-1.22.66-1.75 2-4 .45-.75 1-2.08 1-3 0-2.81-2.48-5-5.5-5zm3.64 7.48c-.25.44-.47.8-.67 1.11-.86 1.41-1.25 2.06-1.45 3.23-.02.05-.02.11-.02.17H5c0-.06 0-.13-.02-.17-.2-1.17-.59-1.83-1.45-3.23-.2-.31-.42-.67-.67-1.11C2.44 6.78 2 5.65 2 5c0-2.2 2.02-4 4.5-4 1.22 0 2.36.42 3.22 1.19C10.55 2.94 11 3.94 11 5c0 .66-.44 1.78-.86 2.48zM4 14h5c-.23 1.14-1.3 2-2.5 2s-2.27-.86-2.5-2z"}))),"Illumina-specific meta information")),(0,l.kt)("div",{parentName:"div",className:"admonition-content"},(0,l.kt)("p",{parentName:"div"},"In addition to the above attributes, there are several attributes specific to Illumina control values that are not described here: Biotin, DNP, Extension, Hyb, NP, NSB, Restore, String, TargetRemoval."))),(0,l.kt)("h2",{id:"header-line-columns"},"Header line columns"),(0,l.kt)("p",null,"Following the meta-information, there are 8 standard header line columns:"),(0,l.kt)("ol",null,(0,l.kt)("li",{parentName:"ol"},"#CHROM"),(0,l.kt)("li",{parentName:"ol"},"POS"),(0,l.kt)("li",{parentName:"ol"},"ID"),(0,l.kt)("li",{parentName:"ol"},"REF"),(0,l.kt)("li",{parentName:"ol"},"ALT"),(0,l.kt)("li",{parentName:"ol"},"QUAL"),(0,l.kt)("li",{parentName:"ol"},"FILTER"),(0,l.kt)("li",{parentName:"ol"},"INFO")),(0,l.kt)("p",null,"These columns are standard to each VCF, but the Illumina Genotyping Array Pipeline has an additional FORMAT column. Pipeline-specific column attributes are described below:"),(0,l.kt)("h3",{id:"filter"},"FILTER"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"DUPE - Filter applied if there are multiple rows in the VCF for the same loci and alleles. For example, if there are two or more rows that share the same chromosome, position, ref allele and alternate alleles, all but one of them will have the \u2018DUPE\u2019 filter set"),(0,l.kt)("li",{parentName:"ul"},"TRIALLELIC - Filter applied if there is a site at which there are two alternate alleles and neither of them is the same as the reference allele"),(0,l.kt)("li",{parentName:"ul"},"ZEROED_OUT_ASSAY - Filter applied if the variant at the site was \u2018zeroed out\u2019 in the Illumina cluster file; typically applied when the calls at the site are intentionally marked as unusual. Genotypes called for sites that are \u2018zeroed out\u2019 will always be no-calls")),(0,l.kt)("h3",{id:"format-genotype"},"FORMAT (genotype)"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"GT - GENOTYPE; standard field described in the ",(0,l.kt)("a",{parentName:"li",href:"https://samtools.github.io/hts-specs/VCFv4.2.pdf"},"VCF specification")),(0,l.kt)("li",{parentName:"ul"},"IGC - Illumina GenCall Confidence Score. Measure of the call confidence"),(0,l.kt)("li",{parentName:"ul"},"X - Raw X intensity as scanned from the original genotyping array"),(0,l.kt)("li",{parentName:"ul"},"Y - Raw Y intensity as scanned from the original genotyping array"),(0,l.kt)("li",{parentName:"ul"},"NORMX - Normalized X intensity"),(0,l.kt)("li",{parentName:"ul"},"NORMY - Normalized Y intensity"),(0,l.kt)("li",{parentName:"ul"},"R - Normalized R Value (one of the polar coordinates after the transformation of NORMX and NORMY)"),(0,l.kt)("li",{parentName:"ul"},"THETA - Normalized Theta value (one of the polar coordinates after the transformation of NORMX and NORMY)"),(0,l.kt)("li",{parentName:"ul"},"LRR - Log R Ratio"),(0,l.kt)("li",{parentName:"ul"},"BAF - B Allele Frequency")),(0,l.kt)("h3",{id:"info"},"INFO"),(0,l.kt)("p",null,"This column contains attributes specific to the array's probe."),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"AC - Allele Count in genotypes for each ALT allele; standard field in the VCF specification"),(0,l.kt)("li",{parentName:"ul"},"AF - Allele Frequency; standard field in the VCF specification"),(0,l.kt)("li",{parentName:"ul"},"AN - Allele Number; standard field in the VCF specification"),(0,l.kt)("li",{parentName:"ul"},"ALLELE_A - Allele A as annotated in the Illumina manifest (a *suffix indicates this is the reference allele)"),(0,l.kt)("li",{parentName:"ul"},"ALLELE_B - Allele B as annotated in the Illumina manifest (a *suffix indicates this is the reference allele)"),(0,l.kt)("li",{parentName:"ul"},"BEADSET_ID - BeadSet ID; an Illumina identifier used for normalization"),(0,l.kt)("li",{parentName:"ul"},"GC_SCORE - Illumina GenTrain Score; a quality score describing the probe design"),(0,l.kt)("li",{parentName:"ul"},"ILLUMINA_BUILD - Genome Build for the design probe sequence, as annotated in the Illumina manifest"),(0,l.kt)("li",{parentName:"ul"},"ILLUMINA_CHR - Chromosome of the design probe sequence as annotated in the Illumina manifest"),(0,l.kt)("li",{parentName:"ul"},"ILLUMINA_POS - Position of the design probe sequence (on ILLUMINA_CHR) as annotated in the Illumina manifest"),(0,l.kt)("li",{parentName:"ul"},"ILLUMINA_STRAND - Strand for the design probe sequence as annotated in the Illumina manifest"),(0,l.kt)("li",{parentName:"ul"},"PROBE_A - Allele A probe sequence as annotated in the Illumina manifest"),(0,l.kt)("li",{parentName:"ul"},"PROBE_B - Allele B probe sequence as annotated in the Illumina manifest; only present on strand ambiguous SNPs"),(0,l.kt)("li",{parentName:"ul"},"SOURCE - Probe source as annotated in the Illumina manifest"),(0,l.kt)("li",{parentName:"ul"},"refSNP - dbSNP rsId for this probe")),(0,l.kt)("p",null,"The remaining attributes describe the cluster definitions provided in the cluster file used for calling the genotype for a SNP."),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"N_AA - Number of AA calls in training set"),(0,l.kt)("li",{parentName:"ul"},"N_AB - Number of AB calls in training set"),(0,l.kt)("li",{parentName:"ul"},"N_BB - Number of BB calls in training set"),(0,l.kt)("li",{parentName:"ul"},"devR_AA - Standard deviation of normalized R for AA cluster"),(0,l.kt)("li",{parentName:"ul"},"devR_AB - Standard deviation of normalized R for AB cluster"),(0,l.kt)("li",{parentName:"ul"},"devR_BB - Standard deviation of normalized R for BB cluster"),(0,l.kt)("li",{parentName:"ul"},"devTHETA_AA - Standard deviation of normalized THETA for AA cluster"),(0,l.kt)("li",{parentName:"ul"},"devTHETA_AB - Standard deviation of normalized THETA for AB cluster"),(0,l.kt)("li",{parentName:"ul"},"devTHETA_BB - Standard deviation of normalized THETA for BB cluster"),(0,l.kt)("li",{parentName:"ul"},"devX_AA - Standard deviation of normalized X for AA cluster"),(0,l.kt)("li",{parentName:"ul"},"devX_AB - Standard deviation of normalized X for AB cluster"),(0,l.kt)("li",{parentName:"ul"},"devX_BB - Standard deviation of normalized X for BB cluster"),(0,l.kt)("li",{parentName:"ul"},"devY_AA - Standard deviation of normalized Y for AA cluster"),(0,l.kt)("li",{parentName:"ul"},"devY_AB - Standard deviation of normalized Y for AB cluster"),(0,l.kt)("li",{parentName:"ul"},"devY_BB - Standard deviation of normalized Y for BB cluster"),(0,l.kt)("li",{parentName:"ul"},"meanR_AA - Mean of normalized R for AA cluster"),(0,l.kt)("li",{parentName:"ul"},"meanR_AB - Mean of normalized R for AB cluster"),(0,l.kt)("li",{parentName:"ul"},"meanR_BB - Mean of normalized R for BB cluster"),(0,l.kt)("li",{parentName:"ul"},"meanTHETA_AA - Mean of normalized THETA for AA cluster"),(0,l.kt)("li",{parentName:"ul"},"meanTHETA_AB - Mean of normalized THETA for AB cluster"),(0,l.kt)("li",{parentName:"ul"},"meanTHETA_BB - Mean of normalized THETA for BB cluster"),(0,l.kt)("li",{parentName:"ul"},"meanX_AA - Mean of normalized X for AA cluster"),(0,l.kt)("li",{parentName:"ul"},"meanX_AB - Mean of normalized X for AB cluster"),(0,l.kt)("li",{parentName:"ul"},"meanX_BB - Mean of normalized X for BB cluster"),(0,l.kt)("li",{parentName:"ul"},"meanY_AA - Mean of normalized Y for AA cluster"),(0,l.kt)("li",{parentName:"ul"},"meanY_AB - Mean of normalized Y for AB cluster"),(0,l.kt)("li",{parentName:"ul"},"meanY_BB - Mean of normalized Y for BB cluster")))}u.isMDXComponent=!0}}]);