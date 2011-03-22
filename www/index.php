
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
    </ul>
  </ul>
</ul>
<hr>
<b>Download:</b>
<blockquote>
  <blockquote>
    <blockquote> 
    </blockquote>
  </blockquote>
</blockquote>
<hr><b>Documentation:</b>
<blockquote>
  <blockquote>
    <ul>
      <li>
    </ul>
  </blockquote>
</blockquote>
<hr><b>Developers:</b>
<blockquote>
  <blockquote>
    <ul>
    </ul>
  </blockquote>
</blockquote>
<hr><b>Financial support:</b>
<blockquote>
  <blockquote> The development of the pacakge was financially supported by the German Federal Ministry of Education and Research (BMBF)
within the AgroClustEr Synbreed .. Synergistic plant and animal breeding (FKZ 0315528A)
    <ul>
    </ul>
  </blockquote>
</blockquote>
<!-- 
<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>
-->
</body>
</html>
