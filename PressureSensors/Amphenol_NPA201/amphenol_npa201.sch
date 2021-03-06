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
<layer number="2" name="Route2" color="1" fill="3" visible="no" active="no"/>
<layer number="3" name="Route3" color="4" fill="3" visible="no" active="no"/>
<layer number="4" name="Route4" color="1" fill="4" visible="no" active="no"/>
<layer number="5" name="Route5" color="4" fill="4" visible="no" active="no"/>
<layer number="6" name="Route6" color="1" fill="8" visible="no" active="no"/>
<layer number="7" name="Route7" color="4" fill="8" visible="no" active="no"/>
<layer number="8" name="Route8" color="1" fill="2" visible="no" active="no"/>
<layer number="9" name="Route9" color="4" fill="2" visible="no" active="no"/>
<layer number="10" name="Route10" color="1" fill="7" visible="no" active="no"/>
<layer number="11" name="Route11" color="4" fill="7" visible="no" active="no"/>
<layer number="12" name="Route12" color="1" fill="5" visible="no" active="no"/>
<layer number="13" name="Route13" color="4" fill="5" visible="no" active="no"/>
<layer number="14" name="Route14" color="1" fill="6" visible="no" active="no"/>
<layer number="15" name="Route15" color="4" fill="6" visible="no" active="no"/>
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
<symbol name="+3V3" library_version="1">
<wire x1="1.27" y1="-1.905" x2="0" y2="0" width="0.254" layer="94"/>
<wire x1="0" y1="0" x2="-1.27" y2="-1.905" width="0.254" layer="94"/>
<text x="-2.54" y="-5.08" size="1.778" layer="96" rot="R90">&gt;VALUE</text>
<pin name="+3V3" x="0" y="-2.54" visible="off" length="short" direction="sup" rot="R90"/>
</symbol>
<symbol name="GND" library_version="1">
<wire x1="-1.905" y1="0" x2="1.905" y2="0" width="0.254" layer="94"/>
<text x="-2.54" y="-2.54" size="1.778" layer="96">&gt;VALUE</text>
<pin name="GND" x="0" y="2.54" visible="off" length="short" direction="sup" rot="R270"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="+3V3" prefix="+3V3" library_version="1">
<description>&lt;b&gt;SUPPLY SYMBOL&lt;/b&gt;</description>
<gates>
<gate name="G$1" symbol="+3V3" x="0" y="0"/>
</gates>
<devices>
<device name="">
<technologies>
<technology name=""/>
</technologies>
</device>
</devices>
</deviceset>
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
<library name="wirepad_08-045">
<packages>
<package name="WIREPAD_0,8/0,45">
<pad name="P$1" x="0" y="0" drill="0.45" diameter="0.7096"/>
</package>
</packages>
<symbols>
<symbol name="WIREPAD_0,8/0,45">
<pin name="P$1" x="-5.08" y="0" length="middle"/>
<wire x1="-7.62" y1="2.54" x2="-2.54" y2="-2.54" width="0.254" layer="94"/>
<wire x1="-7.62" y1="-2.54" x2="-2.54" y2="2.54" width="0.254" layer="94"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="WIREPAD_0,8/0,45">
<gates>
<gate name="G$1" symbol="WIREPAD_0,8/0,45" x="2.54" y="0"/>
</gates>
<devices>
<device name="" package="WIREPAD_0,8/0,45">
<connects>
<connect gate="G$1" pin="P$1" pad="P$1"/>
</connects>
<technologies>
<technology name=""/>
</technologies>
</device>
</devices>
</deviceset>
</devicesets>
</library>
<library name="Amphenol_NPA201">
<packages>
<package name="AMPHENOL_NPA201">
<smd name="P$1" x="0.975" y="-0.725" dx="0.4" dy="0.4" layer="1" rot="R180"/>
<smd name="P$2" x="0.325" y="-0.725" dx="0.4" dy="0.4" layer="1"/>
<smd name="P$3" x="-0.325" y="-0.725" dx="0.4" dy="0.4" layer="1"/>
<smd name="P$4" x="-0.975" y="-0.725" dx="0.4" dy="0.4" layer="1"/>
<smd name="P$5" x="-0.975" y="0.725" dx="0.4" dy="0.4" layer="1"/>
<smd name="P$6" x="-0.325" y="0.725" dx="0.4" dy="0.4" layer="1"/>
<smd name="P$7" x="0.325" y="0.725" dx="0.4" dy="0.4" layer="1"/>
<smd name="P$8" x="0.975" y="0.725" dx="0.4" dy="0.4" layer="1"/>
<circle x="1" y="-0.193" radius="0.0559" width="0.127" layer="21"/>
<wire x1="1.35" y1="-1.09" x2="-1.325" y2="-1.09" width="0.0762" layer="21"/>
<wire x1="-1.325" y1="-1.09" x2="-1.325" y2="1.1" width="0.0762" layer="21"/>
<wire x1="-1.325" y1="1.1" x2="1.35" y2="1.1" width="0.0762" layer="21"/>
<wire x1="1.35" y1="1.1" x2="1.35" y2="-1.09" width="0.0762" layer="21"/>
</package>
</packages>
<symbols>
<symbol name="AMPHENOL_NPA201">
<pin name="1.GND" x="-12.7" y="5.08" length="middle"/>
<pin name="2.NC" x="-12.7" y="0" length="middle"/>
<pin name="3.SDA" x="-12.7" y="-5.08" length="middle"/>
<pin name="4.SCL" x="-12.7" y="-10.16" length="middle"/>
<pin name="5.NC" x="15.24" y="-10.16" length="middle" rot="R180"/>
<pin name="6.VDD" x="15.24" y="-5.08" length="middle" rot="R180"/>
<pin name="7.GND" x="15.24" y="0" length="middle" rot="R180"/>
<pin name="8.VDD" x="15.24" y="5.08" length="middle" rot="R180"/>
<wire x1="-10.16" y1="7.62" x2="-10.16" y2="-12.7" width="0.254" layer="94"/>
<wire x1="-10.16" y1="-12.7" x2="12.7" y2="-12.7" width="0.254" layer="94"/>
<wire x1="12.7" y1="-12.7" x2="12.7" y2="7.62" width="0.254" layer="94"/>
<wire x1="12.7" y1="7.62" x2="-10.16" y2="7.62" width="0.254" layer="94"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="AMPHENOL_NPA201">
<gates>
<gate name="G$1" symbol="AMPHENOL_NPA201" x="0" y="0"/>
</gates>
<devices>
<device name="" package="AMPHENOL_NPA201">
<connects>
<connect gate="G$1" pin="1.GND" pad="P$1"/>
<connect gate="G$1" pin="2.NC" pad="P$2"/>
<connect gate="G$1" pin="3.SDA" pad="P$3"/>
<connect gate="G$1" pin="4.SCL" pad="P$4"/>
<connect gate="G$1" pin="5.NC" pad="P$5"/>
<connect gate="G$1" pin="6.VDD" pad="P$6"/>
<connect gate="G$1" pin="7.GND" pad="P$7"/>
<connect gate="G$1" pin="8.VDD" pad="P$8"/>
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
<part name="+3V1" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="+3V3" device=""/>
<part name="+3V2" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="+3V3" device=""/>
<part name="GND1" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND2" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="U$2" library="wirepad_08-045" deviceset="WIREPAD_0,8/0,45" device=""/>
<part name="U$3" library="wirepad_08-045" deviceset="WIREPAD_0,8/0,45" device=""/>
<part name="U$4" library="wirepad_08-045" deviceset="WIREPAD_0,8/0,45" device=""/>
<part name="U$5" library="wirepad_08-045" deviceset="WIREPAD_0,8/0,45" device=""/>
<part name="U$1" library="Amphenol_NPA201" deviceset="AMPHENOL_NPA201" device=""/>
</parts>
<sheets>
<sheet>
<plain>
</plain>
<instances>
<instance part="+3V1" gate="G$1" x="25.4" y="7.62" rot="R270"/>
<instance part="+3V2" gate="G$1" x="25.4" y="-2.54" rot="R270"/>
<instance part="GND1" gate="1" x="25.4" y="2.54" rot="R90"/>
<instance part="GND2" gate="1" x="-20.32" y="7.62" rot="R270"/>
<instance part="U$2" gate="G$1" x="-20.32" y="-2.54"/>
<instance part="U$3" gate="G$1" x="-20.32" y="-7.62"/>
<instance part="U$4" gate="G$1" x="17.78" y="15.24" rot="R270"/>
<instance part="U$5" gate="G$1" x="35.56" y="2.54" rot="R180"/>
<instance part="U$1" gate="G$1" x="0" y="2.54"/>
</instances>
<busses>
</busses>
<nets>
<net name="GND" class="0">
<segment>
<pinref part="GND2" gate="1" pin="GND"/>
<wire x1="-17.78" y1="7.62" x2="-12.7" y2="7.62" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="1.GND"/>
</segment>
<segment>
<pinref part="GND1" gate="1" pin="GND"/>
<wire x1="15.24" y1="2.54" x2="20.32" y2="2.54" width="0.1524" layer="91"/>
<wire x1="20.32" y1="2.54" x2="22.86" y2="2.54" width="0.1524" layer="91"/>
<wire x1="20.32" y1="2.54" x2="20.32" y2="5.08" width="0.1524" layer="91"/>
<junction x="20.32" y="2.54"/>
<wire x1="20.32" y1="5.08" x2="40.64" y2="5.08" width="0.1524" layer="91"/>
<pinref part="U$5" gate="G$1" pin="P$1"/>
<wire x1="40.64" y1="5.08" x2="40.64" y2="2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="7.GND"/>
</segment>
</net>
<net name="N$1" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$1"/>
<wire x1="-12.7" y1="-2.54" x2="-25.4" y2="-2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="3.SDA"/>
</segment>
</net>
<net name="N$2" class="0">
<segment>
<pinref part="U$3" gate="G$1" pin="P$1"/>
<wire x1="-25.4" y1="-7.62" x2="-12.7" y2="-7.62" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="4.SCL"/>
</segment>
</net>
<net name="+3V3" class="0">
<segment>
<pinref part="+3V2" gate="G$1" pin="+3V3"/>
<wire x1="15.24" y1="-2.54" x2="22.86" y2="-2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="6.VDD"/>
</segment>
<segment>
<pinref part="+3V1" gate="G$1" pin="+3V3"/>
<wire x1="15.24" y1="7.62" x2="17.78" y2="7.62" width="0.1524" layer="91"/>
<pinref part="U$4" gate="G$1" pin="P$1"/>
<wire x1="17.78" y1="7.62" x2="22.86" y2="7.62" width="0.1524" layer="91"/>
<wire x1="17.78" y1="20.32" x2="17.78" y2="7.62" width="0.1524" layer="91"/>
<junction x="17.78" y="7.62"/>
<pinref part="U$1" gate="G$1" pin="8.VDD"/>
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
