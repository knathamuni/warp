"use strict";(self.webpackChunkwebsite_2=self.webpackChunkwebsite_2||[]).push([[6099],{3905:function(e,t,n){n.d(t,{Zo:function(){return c},kt:function(){return f}});var r=n(7294);function i(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function o(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);t&&(r=r.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,r)}return n}function a(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?o(Object(n),!0).forEach((function(t){i(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):o(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function s(e,t){if(null==e)return{};var n,r,i=function(e,t){if(null==e)return{};var n,r,i={},o=Object.keys(e);for(r=0;r<o.length;r++)n=o[r],t.indexOf(n)>=0||(i[n]=e[n]);return i}(e,t);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(r=0;r<o.length;r++)n=o[r],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(i[n]=e[n])}return i}var l=r.createContext({}),p=function(e){var t=r.useContext(l),n=t;return e&&(n="function"==typeof e?e(t):a(a({},t),e)),n},c=function(e){var t=p(e.components);return r.createElement(l.Provider,{value:t},e.children)},d="mdxType",u={inlineCode:"code",wrapper:function(e){var t=e.children;return r.createElement(r.Fragment,{},t)}},m=r.forwardRef((function(e,t){var n=e.components,i=e.mdxType,o=e.originalType,l=e.parentName,c=s(e,["components","mdxType","originalType","parentName"]),d=p(n),m=i,f=d["".concat(l,".").concat(m)]||d[m]||u[m]||o;return n?r.createElement(f,a(a({ref:t},c),{},{components:n})):r.createElement(f,a({ref:t},c))}));function f(e,t){var n=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var o=n.length,a=new Array(o);a[0]=m;var s={};for(var l in t)hasOwnProperty.call(t,l)&&(s[l]=t[l]);s.originalType=e,s[d]="string"==typeof e?e:i,a[1]=s;for(var p=2;p<o;p++)a[p]=n[p];return r.createElement.apply(null,a)}return r.createElement.apply(null,n)}m.displayName="MDXCreateElement"},7287:function(e,t,n){n.r(t),n.d(t,{assets:function(){return l},contentTitle:function(){return a},default:function(){return u},frontMatter:function(){return o},metadata:function(){return s},toc:function(){return p}});var r=n(3117),i=(n(7294),n(3905));const o={sidebar_position:2},a="scATAC v1.2.0 Methods",s={unversionedId:"Pipelines/Single_Cell_ATAC_Seq_Pipeline/scatac.methods",id:"Pipelines/Single_Cell_ATAC_Seq_Pipeline/scatac.methods",title:"scATAC v1.2.0 Methods",description:"Below we provide a sample methods section for a publication. For the complete pipeline documentation, see the scATAC Overview.",source:"@site/docs/Pipelines/Single_Cell_ATAC_Seq_Pipeline/scatac.methods.md",sourceDirName:"Pipelines/Single_Cell_ATAC_Seq_Pipeline",slug:"/Pipelines/Single_Cell_ATAC_Seq_Pipeline/scatac.methods",permalink:"/warp/docs/Pipelines/Single_Cell_ATAC_Seq_Pipeline/scatac.methods",draft:!1,editUrl:"https://github.com/broadinstitute/warp/edit/develop/website/docs/Pipelines/Single_Cell_ATAC_Seq_Pipeline/scatac.methods.md",tags:[],version:"current",lastUpdatedBy:"dependabot[bot]",lastUpdatedAt:1700148680,formattedLastUpdatedAt:"Nov 16, 2023",sidebarPosition:2,frontMatter:{sidebar_position:2},sidebar:"docsSidebar",previous:{title:"scATAC Overview",permalink:"/warp/docs/Pipelines/Single_Cell_ATAC_Seq_Pipeline/README"},next:{title:"Slide-seq Overview",permalink:"/warp/docs/Pipelines/SlideSeq_Pipeline/README"}},l={},p=[{value:"Methods",id:"methods",level:2}],c={toc:p},d="wrapper";function u(e){let{components:t,...n}=e;return(0,i.kt)(d,(0,r.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("h1",{id:"scatac-v120-methods"},"scATAC v1.2.0 Methods"),(0,i.kt)("p",null,"Below we provide a sample methods section for a publication. For the complete pipeline documentation, see the ",(0,i.kt)("a",{parentName:"p",href:"/warp/docs/Pipelines/Single_Cell_ATAC_Seq_Pipeline/README"},"scATAC Overview"),"."),(0,i.kt)("h2",{id:"methods"},"Methods"),(0,i.kt)("p",null,"Preprocessing and cell by bin matrix construction were performed using the scATAC v1.2.0 Pipeline (RRID:SCR_018919)."),(0,i.kt)("p",null,"Prior to processing, paired-end FASTQ files were modified with a custom python script (available at ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/r3fang/SnapTools/blob/master/snaptools/dex_fastq.py"},"https://github.com/r3fang/SnapTools/blob/master/snaptools/dex_fastq.py"),") so that readnames were appended with cell barcodes."),(0,i.kt)("p",null,"The appended reads were then aligned to the Hg38 genomic reference using BWA v0.7.17. The resulting aligned BAM was converted into fragments and filtered using the ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/r3fang/SnapTools"},"SnapTools v1.4.7")," SnapPre function with default parameters."),(0,i.kt)("p",null,"The snap-add-bmat function was then used to add cell-by-bin matrices to the resulting Snap file. 10 kb was selected as the default value for the bin size based on the SnapTools recommendation for mid-size datasets. It can be changed by specifying the desired size as an input to this workflow."),(0,i.kt)("p",null,"Custom python scripts were then used to make a ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/broadinstitute/warp-tools/blob/develop/tools/scripts/makeCompliantBAM.py"},"GA4GH-compliant BAM")," and to export select Snap file metrics to ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/broadinstitute/warp-tools/blob/develop/tools/scripts/breakoutSnap.py"},"individual text files"),"."),(0,i.kt)("p",null,"An example of the pipeline and its outputs is available on ",(0,i.kt)("a",{parentName:"p",href:"https://app.terra.bio/#workspaces/brain-initiative-bcdc/SnapATAC_Pipeline"},"Terra")," and more documentation can be found at ",(0,i.kt)("a",{parentName:"p",href:"/warp/docs/Pipelines/Single_Cell_ATAC_Seq_Pipeline/README"},"here"),".\nExamples of genomic reference files and other inputs can be found in the pipeline\u2019s ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/scATAC/example_inputs/human_example.json"},"example JSON"),"."))}u.isMDXComponent=!0}}]);