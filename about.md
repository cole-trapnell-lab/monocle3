---
layout: default
permalink: /
tags: [Monocle, RNA-Seq]
modified: 2014-11-10
---

Monocle is a toolkit for analyzing single-cell gene expression experiments. Monocle was designed for RNA-Seq, but can also work with single cell qPCR. It performs differential expression analysis, and can find genes that differ between cell types or between cell states. When used to study an ongoing biological process such as cell differentiation, Monocle learns that process and places cells in order according to their progress through it. Monocle finds genes that are dynamically regulated during that process.

Monocle was written by Cole Trapnell with input from Davide Cacchiarelli.

Monocle is provided under the OSI-approved Artistic License (version 2.0)

# News

*To get the latest updates on the Monocle project and the rest of the "Tuxedo tools", please subscribe to our [**mailing list**](https://lists.sourceforge.net/lists/listinfo/bowtie-bio-announce)*

<ul class="post-list">
{% for post in site.posts %}
  <li><article><a href="{{ site.url }}{{ post.url }}">{{ post.title }} <span class="entry-date"><time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %d, %Y" }}</time></span></a></article></li>
{% endfor %}
</ul>
