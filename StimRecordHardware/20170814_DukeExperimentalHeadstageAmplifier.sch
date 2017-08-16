<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE eagle SYSTEM "eagle.dtd">
<eagle version="8.2.2">
<drawing>
<settings>
<setting alwaysvectorfont="no"/>
<setting verticaltext="up"/>
</settings>
<grid distance="0.1" unitdist="inch" unit="inch" style="lines" multiple="1" display="no" altdistance="0.01" altunitdist="inch" altunit="inch"/>
<layers>
<layer number="1" name="Top" color="4" fill="1" visible="no" active="no"/>
<layer number="16" name="Bottom" color="1" fill="1" visible="no" active="no"/>
<layer number="17" name="Pads" color="2" fill="1" visible="no" active="no"/>
<layer number="18" name="Vias" color="2" fill="1" visible="no" active="no"/>
<layer number="19" name="Unrouted" color="6" fill="1" visible="no" active="no"/>
<layer number="20" name="Dimension" color="15" fill="1" visible="no" active="no"/>
<layer number="21" name="tPlace" color="7" fill="1" visible="no" active="no"/>
<layer number="22" name="bPlace" color="7" fill="1" visible="no" active="no"/>
<layer number="23" name="tOrigins" color="15" fill="1" visible="no" active="no"/>
<layer number="24" name="bOrigins" color="15" fill="1" visible="no" active="no"/>
<layer number="25" name="tNames" color="7" fill="1" visible="no" active="no"/>
<layer number="26" name="bNames" color="7" fill="1" visible="no" active="no"/>
<layer number="27" name="tValues" color="7" fill="1" visible="no" active="no"/>
<layer number="28" name="bValues" color="7" fill="1" visible="no" active="no"/>
<layer number="29" name="tStop" color="7" fill="3" visible="no" active="no"/>
<layer number="30" name="bStop" color="7" fill="6" visible="no" active="no"/>
<layer number="31" name="tCream" color="7" fill="4" visible="no" active="no"/>
<layer number="32" name="bCream" color="7" fill="5" visible="no" active="no"/>
<layer number="33" name="tFinish" color="6" fill="3" visible="no" active="no"/>
<layer number="34" name="bFinish" color="6" fill="6" visible="no" active="no"/>
<layer number="35" name="tGlue" color="7" fill="4" visible="no" active="no"/>
<layer number="36" name="bGlue" color="7" fill="5" visible="no" active="no"/>
<layer number="37" name="tTest" color="7" fill="1" visible="no" active="no"/>
<layer number="38" name="bTest" color="7" fill="1" visible="no" active="no"/>
<layer number="39" name="tKeepout" color="4" fill="11" visible="no" active="no"/>
<layer number="40" name="bKeepout" color="1" fill="11" visible="no" active="no"/>
<layer number="41" name="tRestrict" color="4" fill="10" visible="no" active="no"/>
<layer number="42" name="bRestrict" color="1" fill="10" visible="no" active="no"/>
<layer number="43" name="vRestrict" color="2" fill="10" visible="no" active="no"/>
<layer number="44" name="Drills" color="7" fill="1" visible="no" active="no"/>
<layer number="45" name="Holes" color="7" fill="1" visible="no" active="no"/>
<layer number="46" name="Milling" color="3" fill="1" visible="no" active="no"/>
<layer number="47" name="Measures" color="7" fill="1" visible="no" active="no"/>
<layer number="48" name="Document" color="7" fill="1" visible="no" active="no"/>
<layer number="49" name="Reference" color="7" fill="1" visible="no" active="no"/>
<layer number="51" name="tDocu" color="7" fill="1" visible="no" active="no"/>
<layer number="52" name="bDocu" color="7" fill="1" visible="no" active="no"/>
<layer number="90" name="Modules" color="5" fill="1" visible="yes" active="yes"/>
<layer number="91" name="Nets" color="2" fill="1" visible="yes" active="yes"/>
<layer number="92" name="Busses" color="1" fill="1" visible="yes" active="yes"/>
<layer number="93" name="Pins" color="2" fill="1" visible="no" active="yes"/>
<layer number="94" name="Symbols" color="4" fill="1" visible="yes" active="yes"/>
<layer number="95" name="Names" color="7" fill="1" visible="yes" active="yes"/>
<layer number="96" name="Values" color="7" fill="1" visible="yes" active="yes"/>
<layer number="97" name="Info" color="7" fill="1" visible="yes" active="yes"/>
<layer number="98" name="Guide" color="6" fill="1" visible="yes" active="yes"/>
</layers>
<schematic xreflabel="%F%N/%S.%C%R" xrefpart="/%S.%C%R">
<libraries>
<library name="supply1" urn="urn:adsk.eagle:library:371">
<description>&lt;b&gt;Supply Symbols&lt;/b&gt;&lt;p&gt;
 GND, VCC, 0V, +5V, -5V, etc.&lt;p&gt;
 Please keep in mind, that these devices are necessary for the
 automatic wiring of the supply signals.&lt;p&gt;
 The pin name defined in the symbol is identical to the net which is to be wired automatically.&lt;p&gt;
 In this library the device names are the same as the pin names of the symbols, therefore the correct signal names appear next to the supply symbols in the schematic.&lt;p&gt;
 &lt;author&gt;Created by librarian@cadsoft.de&lt;/author&gt;</description>
<packages>
</packages>
<symbols>
<symbol name="GND" library_version="1">
<wire x1="-1.905" y1="0" x2="1.905" y2="0" width="0.254" layer="94"/>
<text x="-2.54" y="-2.54" size="1.778" layer="96">&gt;VALUE</text>
<pin name="GND" x="0" y="2.54" visible="off" length="short" direction="sup" rot="R270"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="GND" prefix="GND" library_version="1">
<description>&lt;b&gt;SUPPLY SYMBOL&lt;/b&gt;</description>
<gates>
<gate name="1" symbol="GND" x="0" y="0"/>
</gates>
<devices>
<device name="">
<technologies>
<technology name=""/>
</technologies>
</device>
</devices>
</deviceset>
</devicesets>
</library>
<library name="MIS_38pin">
<packages>
<package name="MIS-019-01-H-D-EM2">
<smd name="GND3" x="0" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P32" x="0.635" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P30" x="1.27" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P28" x="1.905" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P26" x="2.54" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P24" x="3.175" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P22" x="3.81" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P20" x="4.445" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P18" x="5.08" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="GND2" x="5.715" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P15" x="6.35" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P13" x="6.985" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P11" x="7.62" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P9" x="8.255" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P7" x="8.89" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P5" x="9.525" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P3" x="10.16" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="P1" x="10.795" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="GND1" x="11.43" y="-0.762" dx="1.524" dy="0.3048" layer="1" rot="R90"/>
<smd name="GND4" x="11.43" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P2" x="10.795" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P4" x="10.16" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P6" x="9.525" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P8" x="8.89" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P10" x="8.255" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P12" x="7.62" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P14" x="6.985" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P16" x="6.35" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P17" x="5.715" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P19" x="5.08" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P21" x="4.445" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P23" x="3.81" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P25" x="3.175" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P27" x="2.54" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P29" x="1.905" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="P31" x="1.27" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="GND5" x="0.635" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<smd name="GND6" x="0" y="-0.762" dx="1.524" dy="0.3048" layer="16" rot="R90"/>
<text x="10.4394" y="0.127" size="0.508" layer="21">01</text>
<wire x1="0" y1="-0.7874" x2="11.43" y2="-0.7874" width="0.127" layer="21"/>
<wire x1="11.43" y1="-0.7874" x2="11.938" y2="-1.7018" width="0.127" layer="21"/>
<wire x1="-0.508" y1="-1.7018" x2="0" y2="-0.7874" width="0.127" layer="21"/>
<wire x1="11.938" y1="-1.7018" x2="-0.508" y2="-1.7018" width="0.127" layer="21"/>
</package>
</packages>
<symbols>
<symbol name="MIS-019-01-H-D-EM2">
<pin name="GND1" x="-12.7" y="40.64" length="middle"/>
<pin name="P1" x="-12.7" y="35.56" length="middle"/>
<pin name="P3" x="-12.7" y="30.48" length="middle"/>
<pin name="P5" x="-12.7" y="25.4" length="middle"/>
<pin name="P7" x="-12.7" y="20.32" length="middle"/>
<pin name="P9" x="-12.7" y="15.24" length="middle"/>
<pin name="P11" x="-12.7" y="10.16" length="middle"/>
<pin name="P13" x="-12.7" y="5.08" length="middle"/>
<pin name="P15" x="-12.7" y="0" length="middle"/>
<pin name="GND2" x="-12.7" y="-5.08" length="middle"/>
<pin name="P18" x="-12.7" y="-10.16" length="middle"/>
<pin name="P20" x="-12.7" y="-15.24" length="middle"/>
<pin name="P22" x="-12.7" y="-20.32" length="middle"/>
<pin name="P24" x="-12.7" y="-25.4" length="middle"/>
<pin name="P26" x="-12.7" y="-30.48" length="middle"/>
<pin name="P28" x="-12.7" y="-35.56" length="middle"/>
<pin name="P30" x="-12.7" y="-40.64" length="middle"/>
<pin name="P32" x="-12.7" y="-45.72" length="middle"/>
<pin name="GND3" x="-12.7" y="-50.8" length="middle"/>
<pin name="GND6" x="15.24" y="-50.8" length="middle" rot="R180"/>
<pin name="GND5" x="15.24" y="-45.72" length="middle" rot="R180"/>
<pin name="P31" x="15.24" y="-40.64" length="middle" rot="R180"/>
<pin name="P29" x="15.24" y="-35.56" length="middle" rot="R180"/>
<pin name="P27" x="15.24" y="-30.48" length="middle" rot="R180"/>
<pin name="P25" x="15.24" y="-25.4" length="middle" rot="R180"/>
<pin name="P23" x="15.24" y="-20.32" length="middle" rot="R180"/>
<pin name="P21" x="15.24" y="-15.24" length="middle" rot="R180"/>
<pin name="P19" x="15.24" y="-10.16" length="middle" rot="R180"/>
<pin name="P17" x="15.24" y="-5.08" length="middle" rot="R180"/>
<pin name="P16" x="15.24" y="0" length="middle" rot="R180"/>
<pin name="P14" x="15.24" y="5.08" length="middle" rot="R180"/>
<pin name="P12" x="15.24" y="10.16" length="middle" rot="R180"/>
<pin name="P10" x="15.24" y="15.24" length="middle" rot="R180"/>
<pin name="P8" x="15.24" y="20.32" length="middle" rot="R180"/>
<pin name="P6" x="15.24" y="25.4" length="middle" rot="R180"/>
<pin name="P4" x="15.24" y="30.48" length="middle" rot="R180"/>
<pin name="P2" x="15.24" y="35.56" length="middle" rot="R180"/>
<pin name="GND4" x="15.24" y="40.64" length="middle" rot="R180"/>
<wire x1="-10.16" y1="43.18" x2="-10.16" y2="-53.34" width="0.254" layer="94"/>
<wire x1="-10.16" y1="-53.34" x2="12.7" y2="-53.34" width="0.254" layer="94"/>
<wire x1="12.7" y1="-53.34" x2="12.7" y2="43.18" width="0.254" layer="94"/>
<wire x1="12.7" y1="43.18" x2="-10.16" y2="43.18" width="0.254" layer="94"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="MIS-019-01-H-D-EM2">
<gates>
<gate name="G$1" symbol="MIS-019-01-H-D-EM2" x="0" y="0"/>
</gates>
<devices>
<device name="" package="MIS-019-01-H-D-EM2">
<connects>
<connect gate="G$1" pin="GND1" pad="GND1"/>
<connect gate="G$1" pin="GND2" pad="GND2"/>
<connect gate="G$1" pin="GND3" pad="GND3"/>
<connect gate="G$1" pin="GND4" pad="GND4"/>
<connect gate="G$1" pin="GND5" pad="GND5"/>
<connect gate="G$1" pin="GND6" pad="GND6"/>
<connect gate="G$1" pin="P1" pad="P1"/>
<connect gate="G$1" pin="P10" pad="P10"/>
<connect gate="G$1" pin="P11" pad="P11"/>
<connect gate="G$1" pin="P12" pad="P12"/>
<connect gate="G$1" pin="P13" pad="P13"/>
<connect gate="G$1" pin="P14" pad="P14"/>
<connect gate="G$1" pin="P15" pad="P15"/>
<connect gate="G$1" pin="P16" pad="P16"/>
<connect gate="G$1" pin="P17" pad="P17"/>
<connect gate="G$1" pin="P18" pad="P18"/>
<connect gate="G$1" pin="P19" pad="P19"/>
<connect gate="G$1" pin="P2" pad="P2"/>
<connect gate="G$1" pin="P20" pad="P20"/>
<connect gate="G$1" pin="P21" pad="P21"/>
<connect gate="G$1" pin="P22" pad="P22"/>
<connect gate="G$1" pin="P23" pad="P23"/>
<connect gate="G$1" pin="P24" pad="P24"/>
<connect gate="G$1" pin="P25" pad="P25"/>
<connect gate="G$1" pin="P26" pad="P26"/>
<connect gate="G$1" pin="P27" pad="P27"/>
<connect gate="G$1" pin="P28" pad="P28"/>
<connect gate="G$1" pin="P29" pad="P29"/>
<connect gate="G$1" pin="P3" pad="P3"/>
<connect gate="G$1" pin="P30" pad="P30"/>
<connect gate="G$1" pin="P31" pad="P31"/>
<connect gate="G$1" pin="P32" pad="P32"/>
<connect gate="G$1" pin="P4" pad="P4"/>
<connect gate="G$1" pin="P5" pad="P5"/>
<connect gate="G$1" pin="P6" pad="P6"/>
<connect gate="G$1" pin="P7" pad="P7"/>
<connect gate="G$1" pin="P8" pad="P8"/>
<connect gate="G$1" pin="P9" pad="P9"/>
</connects>
<technologies>
<technology name=""/>
</technologies>
</device>
</devices>
</deviceset>
</devicesets>
</library>
</libraries>
<attributes>
</attributes>
<variantdefs>
</variantdefs>
<classes>
<class number="0" name="default" width="0" drill="0">
</class>
</classes>
<parts>
<part name="GND1" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="U$1" library="MIS_38pin" deviceset="MIS-019-01-H-D-EM2" device=""/>
</parts>
<sheets>
<sheet>
<plain>
</plain>
<instances>
<instance part="GND1" gate="1" x="-40.64" y="71.12" rot="R180"/>
<instance part="U$1" gate="G$1" x="-40.64" y="10.16"/>
</instances>
<busses>
</busses>
<nets>
<net name="GND" class="0">
<segment>
<pinref part="U$1" gate="G$1" pin="GND3"/>
<pinref part="U$1" gate="G$1" pin="GND6"/>
<wire x1="-53.34" y1="-40.64" x2="-38.1" y2="-40.64" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="GND2"/>
<wire x1="-38.1" y1="-40.64" x2="-25.4" y2="-40.64" width="0.1524" layer="91"/>
<wire x1="-53.34" y1="5.08" x2="-38.1" y2="5.08" width="0.1524" layer="91"/>
<wire x1="-38.1" y1="5.08" x2="-38.1" y2="50.8" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="GND4"/>
<wire x1="-38.1" y1="50.8" x2="-25.4" y2="50.8" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="GND1"/>
<wire x1="-38.1" y1="50.8" x2="-53.34" y2="50.8" width="0.1524" layer="91"/>
<junction x="-38.1" y="50.8"/>
<wire x1="-38.1" y1="-40.64" x2="-38.1" y2="5.08" width="0.1524" layer="91"/>
<junction x="-38.1" y="-40.64"/>
<junction x="-38.1" y="5.08"/>
<pinref part="GND1" gate="1" pin="GND"/>
<wire x1="-38.1" y1="50.8" x2="-40.64" y2="50.8" width="0.1524" layer="91"/>
<wire x1="-40.64" y1="50.8" x2="-40.64" y2="68.58" width="0.1524" layer="91"/>
</segment>
</net>
</nets>
</sheet>
</sheets>
</schematic>
</drawing>
<compatibility>
<note version="8.2" severity="warning">
Since Version 8.2, Eagle supports online libraries. The ids
of those online libraries will not be understood (or retained)
with this version.
</note>
</compatibility>
</eagle>
