
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<hr><b>Description:</b>
<ul>
  <ul>
    <ul>
      <li> The R-package <tt>synbreed</tt> provides a framework for the analysis of genomic prediction data (Genomic Selection, GWAS, QTL-mapping) within
an open source software. <br>
      <li> Features </li>  
      <ul>
        <li> Data processing </li>
          <ul>
            <li> Combining raw data sources to a gpData object </li>
            <li>  Conversion from and to class cross in package qtl </li>
            <li> Coding marker data into number of copies of the minor allele </li>
            <li> Preselection of markers according to MAF, % missing values and LD  </li>
            <li> Imputation of missing genotypes by marginal allele distribution, family structure for fully homozygous inbred individuals or flanking markers using Beagle </li>
          </ul>
         <li> Data visualization and analysis </li> 
             <ul>
                <li> Summary method for classes gpData, pedigree and relationshipMatrix </li>
                <li> Marker map representation for low and high density maps </li>
                <li> LD computation as r2 and LD decay visualization as as scatterplot or
stacked histogram  </li>
                <li> Pedigree tree and kinship visualization of relatedness between individuals </li>    
             </ul>
           <li> Statistical models </li>  
             <ul>
                <li> Estimation of pedigree based relationship (additive and dominance) </li>
                <li> Marker based relationship    </li>
                <li> Cross-validation for BLUP, Ridge Regression and Bayesian methods  </li>
             </ul>
            </ul>     
            <li>The package comes with ABSOLUTELY NO WARRANTY; for details
see <a href="http://www.gnu.org/copyleft/gpl.html">http://www.gnu.org/copyleft/gpl.html</a>
(GPL). </li> 
      
    </ul>
  </ul>
</ul>
<hr>
<b>Download:</b>
<blockquote>
  <blockquote>
    <blockquote> 
     <li> Workhshop in G&ouml;ttingen 28/29 March 2011 <br>
         <blockquote> 
         <li> <a href="http://www.computationalgenetics.se/QTLMAS08/QTLMAS/DATA.html"> Data </a>  from XII QTL-MAS Workshop 2008, Uppsala
         <li> <a href="R.zip"> R Code chunks </a> for the workshop </li>
         <li> <a href="slides.pdf"> Slides </a>  from Valentin Wimmer </li>
         <li> <a href="praesentation_Martin_Schlather.pdf"> Slides </a> from Martin Schlather </li>
         <li> <a href="praesentation_RF_for_synbreed_2011_03_22.pdf"> Slides  </a>  from Ulrike Ober  </li>
        </blockquote> 
    </blockquote>
  </blockquote>
</blockquote>
<hr><b>Documentation:</b>
<blockquote>
  <blockquote>
    <ul>
       <li><a href="Poster_synbreed_Rpackage_2011.pdf"> Poster </a> from 11th GABI status seminar, March 15th until 17th 2011 in Potsdam
    </ul>
  </blockquote>
</blockquote>
<hr><b>Developers:</b>
<blockquote>
  <blockquote>
    <ul>  
     <li> <a href="http://wzw.tum.de/plantbreeding/index.php?id=mitarbeiter"> Valentin Wimmer </a>, Chair of Plant Breeding, Technische Universit&auml;t M&uuml;nchen  </li>
     <li> <a href="http://wzw.tum.de/plantbreeding/index.php?id=mitarbeiter"> Theresa Albrecht </a>, Chair of Plant Breeding, Technische Universit&auml;t M&uuml;nchen  </li>
    </ul>
  </blockquote>
</blockquote>
<hr><b>Financial support:</b>
<blockquote>
  <blockquote> The development of the package was financially supported by the German Federal Ministry of Education and Research (BMBF)
within the AgroClustEr ''<a href="http://www.synbreed.tum.de">Synbreed </a> Synergistic plant and animal breeding'' (FKZ 0315528A)
    <ul>
    </ul>
  </blockquote>
</blockquote>
<!-- 
<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>
-->
</body>
</html>
