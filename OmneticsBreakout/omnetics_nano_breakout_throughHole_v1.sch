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
<layer number="150" name="pDim" color="6" fill="1" visible="yes" active="yes"/>
<layer number="151" name="BacRum" color="10" fill="1" visible="yes" active="yes"/>
<layer number="152" name="pCourtYard" color="10" fill="1" visible="yes" active="yes"/>
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
<symbol name="AGND" library_version="1">
<wire x1="-1.905" y1="0" x2="1.905" y2="0" width="0.254" layer="94"/>
<wire x1="-1.0922" y1="-0.508" x2="1.0922" y2="-0.508" width="0.254" layer="94"/>
<text x="-2.54" y="-5.08" size="1.778" layer="96" rot="R90">&gt;VALUE</text>
<pin name="AGND" x="0" y="2.54" visible="off" length="short" direction="sup" rot="R270"/>
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
<deviceset name="AGND" prefix="AGND" library_version="1">
<description>&lt;b&gt;SUPPLY SYMBOL&lt;/b&gt;</description>
<gates>
<gate name="VR1" symbol="AGND" x="0" y="0"/>
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
<library name="Sombeck_lib">
<packages>
<package name="FSTH-118-01-F-D">
<pad name="P$1" x="0" y="0" drill="0.6604"/>
<pad name="P$2" x="1.27" y="0" drill="0.6604"/>
<pad name="P$3" x="2.54" y="0" drill="0.6604"/>
<pad name="P$4" x="3.81" y="0" drill="0.6604"/>
<pad name="P$5" x="5.08" y="0" drill="0.6604"/>
<pad name="P$6" x="6.35" y="0" drill="0.6604"/>
<pad name="P$7" x="7.62" y="0" drill="0.6604"/>
<pad name="P$8" x="8.89" y="0" drill="0.6604"/>
<pad name="P$9" x="10.16" y="0" drill="0.6604"/>
<pad name="P$10" x="11.43" y="0" drill="0.6604"/>
<pad name="P$11" x="12.7" y="0" drill="0.6604"/>
<pad name="P$12" x="13.97" y="0" drill="0.6604"/>
<pad name="P$13" x="15.24" y="0" drill="0.6604"/>
<pad name="P$14" x="16.51" y="0" drill="0.6604"/>
<pad name="P$15" x="17.78" y="0" drill="0.6604"/>
<pad name="P$16" x="19.05" y="0" drill="0.6604"/>
<pad name="P$17" x="20.32" y="0" drill="0.6604"/>
<pad name="P$18" x="21.59" y="0" drill="0.6604"/>
<pad name="P$19" x="21.59" y="1.27" drill="0.6604"/>
<pad name="P$20" x="20.32" y="1.27" drill="0.6604"/>
<pad name="P$21" x="19.05" y="1.27" drill="0.6604"/>
<pad name="P$22" x="17.78" y="1.27" drill="0.6604"/>
<pad name="P$23" x="16.51" y="1.27" drill="0.6604"/>
<pad name="P$24" x="15.24" y="1.27" drill="0.6604"/>
<pad name="P$25" x="13.97" y="1.27" drill="0.6604"/>
<pad name="P$26" x="12.7" y="1.27" drill="0.6604"/>
<pad name="P$27" x="11.43" y="1.27" drill="0.6604"/>
<pad name="P$28" x="10.16" y="1.27" drill="0.6604"/>
<pad name="P$29" x="8.89" y="1.27" drill="0.6604"/>
<pad name="P$30" x="7.62" y="1.27" drill="0.6604"/>
<pad name="P$31" x="6.35" y="1.27" drill="0.6604"/>
<pad name="P$32" x="5.08" y="1.27" drill="0.6604"/>
<pad name="P$33" x="3.81" y="1.27" drill="0.6604"/>
<pad name="P$34" x="2.54" y="1.27" drill="0.6604"/>
<pad name="P$35" x="1.27" y="1.27" drill="0.6604"/>
<pad name="P$36" x="0" y="1.27" drill="0.6604"/>
<wire x1="-1.27" y1="2.54" x2="-1.27" y2="-1.27" width="0.127" layer="21"/>
<wire x1="-1.27" y1="-1.27" x2="22.86" y2="-1.27" width="0.127" layer="21"/>
<wire x1="22.86" y1="-1.27" x2="22.86" y2="2.54" width="0.127" layer="21"/>
<wire x1="22.86" y1="2.54" x2="-1.27" y2="2.54" width="0.127" layer="21"/>
<text x="-2.54" y="-1.27" size="1.27" layer="21">1</text>
</package>
<package name="OMNETICS_NANO_18X2">
<smd name="REF0" x="0" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="GND0" x="0" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_1" x="0.635" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_3" x="1.27" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_5" x="1.905" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_7" x="2.54" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_9" x="3.175" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_11" x="3.81" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_13" x="4.445" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_15" x="5.08" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_2" x="0.635" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_4" x="1.27" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_6" x="1.905" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_8" x="2.54" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_10" x="3.175" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_12" x="3.81" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_14" x="4.445" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_16" x="5.08" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<text x="-0.508" y="-0.3556" size="0.4064" layer="21">1</text>
<smd name="P_18" x="5.715" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_20" x="6.35" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_22" x="6.985" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_24" x="7.62" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_26" x="8.255" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_28" x="8.89" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_30" x="9.525" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_32" x="10.16" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="GND34" x="10.795" y="0" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_17" x="5.715" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_19" x="6.35" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_21" x="6.985" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_23" x="7.62" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_25" x="8.255" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_27" x="8.89" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_29" x="9.525" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="P_31" x="10.16" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
<smd name="REF33" x="10.795" y="1.27" dx="1.016" dy="0.381" layer="1" rot="R90"/>
</package>
<package name="OMNETICS-NANO-THROUGHHOLE">
<pad name="P$1" x="0" y="0.762" drill="0.254" diameter="0.508"/>
<pad name="P$2" x="0" y="0" drill="0.3302" diameter="0.4572"/>
<pad name="P$3" x="0.635" y="0.762" drill="0.254" diameter="0.508"/>
<pad name="P$4" x="0.635" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$5" x="1.27" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$6" x="1.27" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$7" x="1.905" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$8" x="1.905" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$9" x="2.54" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$10" x="2.54" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$11" x="3.175" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$12" x="3.175" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$13" x="3.81" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$14" x="3.81" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$15" x="4.445" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$16" x="4.445" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$17" x="5.08" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$18" x="5.08" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$19" x="5.715" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$20" x="5.715" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$21" x="6.35" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$22" x="6.35" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$23" x="6.985" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$24" x="6.985" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$25" x="7.62" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$26" x="7.62" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$27" x="8.255" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$28" x="8.255" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$29" x="8.89" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$30" x="8.89" y="0" drill="0.35" diameter="0.508"/>
<pad name="P$31" x="9.525" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="P$32" x="9.525" y="0" drill="0.35" diameter="0.508"/>
<pad name="REF33" x="10.16" y="0.762" drill="0.35" diameter="0.508"/>
<pad name="GND34" x="10.16" y="0" drill="0.35" diameter="0.508"/>
<pad name="GND0" x="-0.635" y="0.762" drill="0.254" diameter="0.508"/>
<pad name="REF0" x="-0.635" y="0" drill="0.3302" diameter="0.4572"/>
<wire x1="-1.27" y1="1.27" x2="-1.27" y2="-0.635" width="0.127" layer="21"/>
<wire x1="-1.27" y1="-0.635" x2="10.795" y2="-0.635" width="0.127" layer="21"/>
<wire x1="10.795" y1="-0.635" x2="10.795" y2="1.27" width="0.127" layer="21"/>
<wire x1="10.795" y1="1.27" x2="-1.27" y2="1.27" width="0.127" layer="21"/>
<text x="-1.905" y="0.635" size="0.635" layer="21">1</text>
</package>
</packages>
<symbols>
<symbol name="FSTH-118-01-F-D">
<pin name="P$1" x="-10.16" y="15.24" length="middle"/>
<pin name="P$2" x="-10.16" y="12.7" length="middle"/>
<pin name="P$3" x="-10.16" y="10.16" length="middle"/>
<pin name="P$4" x="-10.16" y="7.62" length="middle"/>
<pin name="P$5" x="-10.16" y="5.08" length="middle"/>
<pin name="P$6" x="-10.16" y="2.54" length="middle"/>
<pin name="P$7" x="-10.16" y="0" length="middle"/>
<pin name="P$8" x="-10.16" y="-2.54" length="middle"/>
<pin name="P$9" x="-10.16" y="-5.08" length="middle"/>
<pin name="P$10" x="-10.16" y="-7.62" length="middle"/>
<pin name="P$11" x="-10.16" y="-10.16" length="middle"/>
<pin name="P$12" x="-10.16" y="-12.7" length="middle"/>
<pin name="P$13" x="-10.16" y="-15.24" length="middle"/>
<pin name="P$14" x="-10.16" y="-17.78" length="middle"/>
<pin name="P$15" x="-10.16" y="-20.32" length="middle"/>
<pin name="P$16" x="-10.16" y="-22.86" length="middle"/>
<pin name="P$17" x="-10.16" y="-25.4" length="middle"/>
<pin name="P$18" x="-10.16" y="-27.94" length="middle"/>
<pin name="P$19" x="15.24" y="-27.94" length="middle" rot="R180"/>
<pin name="P$20" x="15.24" y="-25.4" length="middle" rot="R180"/>
<pin name="P$21" x="15.24" y="-22.86" length="middle" rot="R180"/>
<pin name="P$22" x="15.24" y="-20.32" length="middle" rot="R180"/>
<pin name="P$23" x="15.24" y="-17.78" length="middle" rot="R180"/>
<pin name="P$24" x="15.24" y="-15.24" length="middle" rot="R180"/>
<pin name="P$25" x="15.24" y="-12.7" length="middle" rot="R180"/>
<pin name="P$26" x="15.24" y="-10.16" length="middle" rot="R180"/>
<pin name="P$27" x="15.24" y="-7.62" length="middle" rot="R180"/>
<pin name="P$28" x="15.24" y="-5.08" length="middle" rot="R180"/>
<pin name="P$29" x="15.24" y="-2.54" length="middle" rot="R180"/>
<pin name="P$30" x="15.24" y="0" length="middle" rot="R180"/>
<pin name="P$31" x="15.24" y="2.54" length="middle" rot="R180"/>
<pin name="P$32" x="15.24" y="5.08" length="middle" rot="R180"/>
<pin name="P$33" x="15.24" y="7.62" length="middle" rot="R180"/>
<pin name="P$34" x="15.24" y="10.16" length="middle" rot="R180"/>
<pin name="P$35" x="15.24" y="12.7" length="middle" rot="R180"/>
<pin name="P$36" x="15.24" y="15.24" length="middle" rot="R180"/>
<wire x1="-7.62" y1="17.78" x2="12.7" y2="17.78" width="0.254" layer="94"/>
<wire x1="12.7" y1="17.78" x2="12.7" y2="-30.48" width="0.254" layer="94"/>
<wire x1="12.7" y1="-30.48" x2="-7.62" y2="-30.48" width="0.254" layer="94"/>
<wire x1="-7.62" y1="-30.48" x2="-7.62" y2="17.78" width="0.254" layer="94"/>
</symbol>
<symbol name="OMNETICS_NANO_18X2">
<pin name="GND0" x="-15.24" y="-43.18" length="middle"/>
<pin name="P$1" x="-15.24" y="-38.1" length="middle"/>
<pin name="P$2" x="12.7" y="-38.1" length="middle" rot="R180"/>
<pin name="P$3" x="-15.24" y="-33.02" length="middle"/>
<pin name="P$4" x="12.7" y="-33.02" length="middle" rot="R180"/>
<pin name="P$5" x="-15.24" y="-27.94" length="middle"/>
<pin name="P$6" x="12.7" y="-27.94" length="middle" rot="R180"/>
<pin name="P$7" x="-15.24" y="-22.86" length="middle"/>
<pin name="P$8" x="12.7" y="-22.86" length="middle" rot="R180"/>
<pin name="P$9" x="-15.24" y="-17.78" length="middle"/>
<pin name="P$10" x="12.7" y="-17.78" length="middle" rot="R180"/>
<pin name="P$11" x="-15.24" y="-12.7" length="middle"/>
<pin name="P$12" x="12.7" y="-12.7" length="middle" rot="R180"/>
<pin name="P$13" x="-15.24" y="-7.62" length="middle"/>
<pin name="P$14" x="12.7" y="-7.62" length="middle" rot="R180"/>
<pin name="P$15" x="-15.24" y="-2.54" length="middle"/>
<pin name="P$16" x="12.7" y="-2.54" length="middle" rot="R180"/>
<pin name="P$17" x="-15.24" y="2.54" length="middle"/>
<pin name="P$18" x="12.7" y="2.54" length="middle" rot="R180"/>
<pin name="P$19" x="-15.24" y="7.62" length="middle"/>
<pin name="P$20" x="12.7" y="7.62" length="middle" rot="R180"/>
<pin name="P$21" x="-15.24" y="12.7" length="middle"/>
<pin name="P$22" x="12.7" y="12.7" length="middle" rot="R180"/>
<pin name="P$23" x="-15.24" y="17.78" length="middle"/>
<pin name="P$24" x="12.7" y="17.78" length="middle" rot="R180"/>
<pin name="P$25" x="-15.24" y="22.86" length="middle"/>
<pin name="P$26" x="12.7" y="22.86" length="middle" rot="R180"/>
<pin name="P$27" x="-15.24" y="27.94" length="middle"/>
<pin name="P$28" x="12.7" y="27.94" length="middle" rot="R180"/>
<pin name="P$29" x="-15.24" y="33.02" length="middle"/>
<pin name="P$30" x="12.7" y="33.02" length="middle" rot="R180"/>
<pin name="P$31" x="-15.24" y="38.1" length="middle"/>
<pin name="P$32" x="12.7" y="38.1" length="middle" rot="R180"/>
<pin name="REF0" x="12.7" y="-43.18" length="middle" rot="R180"/>
<pin name="REF33" x="-15.24" y="43.18" length="middle"/>
<pin name="GND34" x="12.7" y="43.18" length="middle" rot="R180"/>
<wire x1="-12.7" y1="45.72" x2="10.16" y2="45.72" width="0.254" layer="94"/>
<wire x1="10.16" y1="45.72" x2="10.16" y2="-45.72" width="0.254" layer="94"/>
<wire x1="10.16" y1="-45.72" x2="-12.7" y2="-45.72" width="0.254" layer="94"/>
<wire x1="-12.7" y1="-45.72" x2="-12.7" y2="45.72" width="0.254" layer="94"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="FSTH-118-01-F-D">
<gates>
<gate name="G$1" symbol="FSTH-118-01-F-D" x="-2.54" y="5.08"/>
</gates>
<devices>
<device name="" package="FSTH-118-01-F-D">
<connects>
<connect gate="G$1" pin="P$1" pad="P$1"/>
<connect gate="G$1" pin="P$10" pad="P$10"/>
<connect gate="G$1" pin="P$11" pad="P$11"/>
<connect gate="G$1" pin="P$12" pad="P$12"/>
<connect gate="G$1" pin="P$13" pad="P$13"/>
<connect gate="G$1" pin="P$14" pad="P$14"/>
<connect gate="G$1" pin="P$15" pad="P$15"/>
<connect gate="G$1" pin="P$16" pad="P$16"/>
<connect gate="G$1" pin="P$17" pad="P$17"/>
<connect gate="G$1" pin="P$18" pad="P$18"/>
<connect gate="G$1" pin="P$19" pad="P$19"/>
<connect gate="G$1" pin="P$2" pad="P$2"/>
<connect gate="G$1" pin="P$20" pad="P$20"/>
<connect gate="G$1" pin="P$21" pad="P$21"/>
<connect gate="G$1" pin="P$22" pad="P$22"/>
<connect gate="G$1" pin="P$23" pad="P$23"/>
<connect gate="G$1" pin="P$24" pad="P$24"/>
<connect gate="G$1" pin="P$25" pad="P$25"/>
<connect gate="G$1" pin="P$26" pad="P$26"/>
<connect gate="G$1" pin="P$27" pad="P$27"/>
<connect gate="G$1" pin="P$28" pad="P$28"/>
<connect gate="G$1" pin="P$29" pad="P$29"/>
<connect gate="G$1" pin="P$3" pad="P$3"/>
<connect gate="G$1" pin="P$30" pad="P$30"/>
<connect gate="G$1" pin="P$31" pad="P$31"/>
<connect gate="G$1" pin="P$32" pad="P$32"/>
<connect gate="G$1" pin="P$33" pad="P$33"/>
<connect gate="G$1" pin="P$34" pad="P$34"/>
<connect gate="G$1" pin="P$35" pad="P$35"/>
<connect gate="G$1" pin="P$36" pad="P$36"/>
<connect gate="G$1" pin="P$4" pad="P$4"/>
<connect gate="G$1" pin="P$5" pad="P$5"/>
<connect gate="G$1" pin="P$6" pad="P$6"/>
<connect gate="G$1" pin="P$7" pad="P$7"/>
<connect gate="G$1" pin="P$8" pad="P$8"/>
<connect gate="G$1" pin="P$9" pad="P$9"/>
</connects>
<technologies>
<technology name=""/>
</technologies>
</device>
</devices>
</deviceset>
<deviceset name="OMNETICS_NANO_18X2">
<gates>
<gate name="G$1" symbol="OMNETICS_NANO_18X2" x="0" y="0"/>
</gates>
<devices>
<device name="LCONNECTOR" package="OMNETICS_NANO_18X2">
<connects>
<connect gate="G$1" pin="GND0" pad="GND0"/>
<connect gate="G$1" pin="GND34" pad="GND34"/>
<connect gate="G$1" pin="P$1" pad="P_1"/>
<connect gate="G$1" pin="P$10" pad="P_10"/>
<connect gate="G$1" pin="P$11" pad="P_11"/>
<connect gate="G$1" pin="P$12" pad="P_12"/>
<connect gate="G$1" pin="P$13" pad="P_13"/>
<connect gate="G$1" pin="P$14" pad="P_14"/>
<connect gate="G$1" pin="P$15" pad="P_15"/>
<connect gate="G$1" pin="P$16" pad="P_16"/>
<connect gate="G$1" pin="P$17" pad="P_17"/>
<connect gate="G$1" pin="P$18" pad="P_18"/>
<connect gate="G$1" pin="P$19" pad="P_19"/>
<connect gate="G$1" pin="P$2" pad="P_2"/>
<connect gate="G$1" pin="P$20" pad="P_20"/>
<connect gate="G$1" pin="P$21" pad="P_21"/>
<connect gate="G$1" pin="P$22" pad="P_22"/>
<connect gate="G$1" pin="P$23" pad="P_23"/>
<connect gate="G$1" pin="P$24" pad="P_24"/>
<connect gate="G$1" pin="P$25" pad="P_25"/>
<connect gate="G$1" pin="P$26" pad="P_26"/>
<connect gate="G$1" pin="P$27" pad="P_27"/>
<connect gate="G$1" pin="P$28" pad="P_28"/>
<connect gate="G$1" pin="P$29" pad="P_29"/>
<connect gate="G$1" pin="P$3" pad="P_3"/>
<connect gate="G$1" pin="P$30" pad="P_30"/>
<connect gate="G$1" pin="P$31" pad="P_31"/>
<connect gate="G$1" pin="P$32" pad="P_32"/>
<connect gate="G$1" pin="P$4" pad="P_4"/>
<connect gate="G$1" pin="P$5" pad="P_5"/>
<connect gate="G$1" pin="P$6" pad="P_6"/>
<connect gate="G$1" pin="P$7" pad="P_7"/>
<connect gate="G$1" pin="P$8" pad="P_8"/>
<connect gate="G$1" pin="P$9" pad="P_9"/>
<connect gate="G$1" pin="REF0" pad="REF0"/>
<connect gate="G$1" pin="REF33" pad="REF33"/>
</connects>
<technologies>
<technology name=""/>
</technologies>
</device>
<device name="" package="OMNETICS-NANO-THROUGHHOLE">
<connects>
<connect gate="G$1" pin="GND0" pad="GND0"/>
<connect gate="G$1" pin="GND34" pad="GND34"/>
<connect gate="G$1" pin="P$1" pad="P$1"/>
<connect gate="G$1" pin="P$10" pad="P$10"/>
<connect gate="G$1" pin="P$11" pad="P$11"/>
<connect gate="G$1" pin="P$12" pad="P$12"/>
<connect gate="G$1" pin="P$13" pad="P$13"/>
<connect gate="G$1" pin="P$14" pad="P$14"/>
<connect gate="G$1" pin="P$15" pad="P$15"/>
<connect gate="G$1" pin="P$16" pad="P$16"/>
<connect gate="G$1" pin="P$17" pad="P$17"/>
<connect gate="G$1" pin="P$18" pad="P$18"/>
<connect gate="G$1" pin="P$19" pad="P$19"/>
<connect gate="G$1" pin="P$2" pad="P$2"/>
<connect gate="G$1" pin="P$20" pad="P$20"/>
<connect gate="G$1" pin="P$21" pad="P$21"/>
<connect gate="G$1" pin="P$22" pad="P$22"/>
<connect gate="G$1" pin="P$23" pad="P$23"/>
<connect gate="G$1" pin="P$24" pad="P$24"/>
<connect gate="G$1" pin="P$25" pad="P$25"/>
<connect gate="G$1" pin="P$26" pad="P$26"/>
<connect gate="G$1" pin="P$27" pad="P$27"/>
<connect gate="G$1" pin="P$28" pad="P$28"/>
<connect gate="G$1" pin="P$29" pad="P$29"/>
<connect gate="G$1" pin="P$3" pad="P$3"/>
<connect gate="G$1" pin="P$30" pad="P$30"/>
<connect gate="G$1" pin="P$31" pad="P$31"/>
<connect gate="G$1" pin="P$32" pad="P$32"/>
<connect gate="G$1" pin="P$4" pad="P$4"/>
<connect gate="G$1" pin="P$5" pad="P$5"/>
<connect gate="G$1" pin="P$6" pad="P$6"/>
<connect gate="G$1" pin="P$7" pad="P$7"/>
<connect gate="G$1" pin="P$8" pad="P$8"/>
<connect gate="G$1" pin="P$9" pad="P$9"/>
<connect gate="G$1" pin="REF0" pad="REF0"/>
<connect gate="G$1" pin="REF33" pad="REF33"/>
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
<part name="GND2" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND3" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND4" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="AGND1" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="AGND" device=""/>
<part name="AGND2" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="AGND" device=""/>
<part name="U$1" library="Sombeck_lib" deviceset="FSTH-118-01-F-D" device=""/>
<part name="U$2" library="Sombeck_lib" deviceset="OMNETICS_NANO_18X2" device=""/>
</parts>
<sheets>
<sheet>
<plain>
</plain>
<instances>
<instance part="GND1" gate="1" x="63.5" y="40.64" rot="R90"/>
<instance part="GND2" gate="1" x="15.24" y="38.1" rot="R270"/>
<instance part="GND3" gate="1" x="-12.7" y="40.64" rot="R90"/>
<instance part="GND4" gate="1" x="-63.5" y="-45.72" rot="R270"/>
<instance part="AGND1" gate="VR1" x="-66.04" y="40.64" rot="R270"/>
<instance part="AGND2" gate="VR1" x="60.96" y="38.1" rot="R90"/>
<instance part="U$1" gate="G$1" x="43.18" y="12.7" rot="R180"/>
<instance part="U$2" gate="G$1" x="-38.1" y="-2.54"/>
</instances>
<busses>
</busses>
<nets>
<net name="GND" class="0">
<segment>
<pinref part="GND1" gate="1" pin="GND"/>
<wire x1="53.34" y1="40.64" x2="60.96" y2="40.64" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$18"/>
</segment>
<segment>
<pinref part="GND2" gate="1" pin="GND"/>
<pinref part="U$1" gate="G$1" pin="P$20"/>
<wire x1="17.78" y1="38.1" x2="27.94" y2="38.1" width="0.1524" layer="91"/>
</segment>
<segment>
<pinref part="GND4" gate="1" pin="GND"/>
<wire x1="-53.34" y1="-45.72" x2="-60.96" y2="-45.72" width="0.1524" layer="91"/>
<pinref part="U$2" gate="G$1" pin="GND0"/>
</segment>
<segment>
<pinref part="GND3" gate="1" pin="GND"/>
<wire x1="-25.4" y1="40.64" x2="-15.24" y2="40.64" width="0.1524" layer="91"/>
<pinref part="U$2" gate="G$1" pin="GND34"/>
</segment>
</net>
<net name="AGND" class="0">
<segment>
<pinref part="AGND1" gate="VR1" pin="AGND"/>
<wire x1="-53.34" y1="40.64" x2="-63.5" y2="40.64" width="0.1524" layer="91"/>
<pinref part="U$2" gate="G$1" pin="REF33"/>
</segment>
<segment>
<pinref part="AGND2" gate="VR1" pin="AGND"/>
<wire x1="58.42" y1="38.1" x2="53.34" y2="38.1" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$17"/>
</segment>
</net>
<net name="N$1" class="0">
<segment>
<wire x1="-53.34" y1="-40.64" x2="-68.58" y2="-40.64" width="0.1524" layer="91"/>
<wire x1="-68.58" y1="-40.64" x2="-68.58" y2="-53.34" width="0.1524" layer="91"/>
<wire x1="-68.58" y1="-53.34" x2="60.96" y2="-53.34" width="0.1524" layer="91"/>
<wire x1="60.96" y1="-53.34" x2="60.96" y2="-2.54" width="0.1524" layer="91"/>
<wire x1="60.96" y1="-2.54" x2="53.34" y2="-2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$1"/>
<pinref part="U$2" gate="G$1" pin="P$1"/>
</segment>
</net>
<net name="N$2" class="0">
<segment>
<wire x1="27.94" y1="-2.54" x2="15.24" y2="-2.54" width="0.1524" layer="91"/>
<wire x1="15.24" y1="-2.54" x2="15.24" y2="-40.64" width="0.1524" layer="91"/>
<wire x1="15.24" y1="-40.64" x2="-25.4" y2="-40.64" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$36"/>
<pinref part="U$2" gate="G$1" pin="P$2"/>
</segment>
</net>
<net name="N$3" class="0">
<segment>
<wire x1="27.94" y1="0" x2="12.7" y2="0" width="0.1524" layer="91"/>
<wire x1="12.7" y1="0" x2="12.7" y2="-35.56" width="0.1524" layer="91"/>
<wire x1="12.7" y1="-35.56" x2="-25.4" y2="-35.56" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$35"/>
<pinref part="U$2" gate="G$1" pin="P$4"/>
</segment>
</net>
<net name="N$4" class="0">
<segment>
<wire x1="-25.4" y1="-30.48" x2="10.16" y2="-30.48" width="0.1524" layer="91"/>
<wire x1="10.16" y1="-30.48" x2="10.16" y2="2.54" width="0.1524" layer="91"/>
<wire x1="10.16" y1="2.54" x2="27.94" y2="2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$34"/>
<pinref part="U$2" gate="G$1" pin="P$6"/>
</segment>
</net>
<net name="N$5" class="0">
<segment>
<wire x1="27.94" y1="5.08" x2="7.62" y2="5.08" width="0.1524" layer="91"/>
<wire x1="7.62" y1="5.08" x2="7.62" y2="-25.4" width="0.1524" layer="91"/>
<wire x1="7.62" y1="-25.4" x2="-25.4" y2="-25.4" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$33"/>
<pinref part="U$2" gate="G$1" pin="P$8"/>
</segment>
</net>
<net name="N$6" class="0">
<segment>
<wire x1="-25.4" y1="-20.32" x2="5.08" y2="-20.32" width="0.1524" layer="91"/>
<wire x1="5.08" y1="-20.32" x2="5.08" y2="7.62" width="0.1524" layer="91"/>
<wire x1="5.08" y1="7.62" x2="27.94" y2="7.62" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$32"/>
<pinref part="U$2" gate="G$1" pin="P$10"/>
</segment>
</net>
<net name="N$7" class="0">
<segment>
<wire x1="27.94" y1="10.16" x2="2.54" y2="10.16" width="0.1524" layer="91"/>
<wire x1="2.54" y1="10.16" x2="2.54" y2="-15.24" width="0.1524" layer="91"/>
<wire x1="2.54" y1="-15.24" x2="-25.4" y2="-15.24" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$31"/>
<pinref part="U$2" gate="G$1" pin="P$12"/>
</segment>
</net>
<net name="N$8" class="0">
<segment>
<wire x1="-25.4" y1="-10.16" x2="0" y2="-10.16" width="0.1524" layer="91"/>
<wire x1="0" y1="-10.16" x2="0" y2="12.7" width="0.1524" layer="91"/>
<wire x1="0" y1="12.7" x2="27.94" y2="12.7" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$30"/>
<pinref part="U$2" gate="G$1" pin="P$14"/>
</segment>
</net>
<net name="N$9" class="0">
<segment>
<wire x1="-25.4" y1="-5.08" x2="-2.54" y2="-5.08" width="0.1524" layer="91"/>
<wire x1="-2.54" y1="-5.08" x2="-2.54" y2="15.24" width="0.1524" layer="91"/>
<wire x1="-2.54" y1="15.24" x2="27.94" y2="15.24" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$29"/>
<pinref part="U$2" gate="G$1" pin="P$16"/>
</segment>
</net>
<net name="N$10" class="0">
<segment>
<wire x1="-25.4" y1="0" x2="-5.08" y2="0" width="0.1524" layer="91"/>
<wire x1="-5.08" y1="0" x2="-5.08" y2="17.78" width="0.1524" layer="91"/>
<wire x1="-5.08" y1="17.78" x2="27.94" y2="17.78" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$28"/>
<pinref part="U$2" gate="G$1" pin="P$18"/>
</segment>
</net>
<net name="N$11" class="0">
<segment>
<wire x1="-25.4" y1="5.08" x2="-7.62" y2="5.08" width="0.1524" layer="91"/>
<wire x1="-7.62" y1="5.08" x2="-7.62" y2="20.32" width="0.1524" layer="91"/>
<wire x1="-7.62" y1="20.32" x2="27.94" y2="20.32" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$27"/>
<pinref part="U$2" gate="G$1" pin="P$20"/>
</segment>
</net>
<net name="N$12" class="0">
<segment>
<wire x1="-25.4" y1="10.16" x2="-10.16" y2="10.16" width="0.1524" layer="91"/>
<wire x1="-10.16" y1="10.16" x2="-10.16" y2="22.86" width="0.1524" layer="91"/>
<wire x1="-10.16" y1="22.86" x2="27.94" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$26"/>
<pinref part="U$2" gate="G$1" pin="P$22"/>
</segment>
</net>
<net name="N$13" class="0">
<segment>
<wire x1="-25.4" y1="15.24" x2="-12.7" y2="15.24" width="0.1524" layer="91"/>
<wire x1="-12.7" y1="15.24" x2="-12.7" y2="25.4" width="0.1524" layer="91"/>
<wire x1="-12.7" y1="25.4" x2="27.94" y2="25.4" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$25"/>
<pinref part="U$2" gate="G$1" pin="P$24"/>
</segment>
</net>
<net name="N$14" class="0">
<segment>
<wire x1="-25.4" y1="20.32" x2="-15.24" y2="20.32" width="0.1524" layer="91"/>
<wire x1="-15.24" y1="20.32" x2="-15.24" y2="27.94" width="0.1524" layer="91"/>
<wire x1="-15.24" y1="27.94" x2="27.94" y2="27.94" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$24"/>
<pinref part="U$2" gate="G$1" pin="P$26"/>
</segment>
</net>
<net name="N$15" class="0">
<segment>
<wire x1="-25.4" y1="25.4" x2="-17.78" y2="25.4" width="0.1524" layer="91"/>
<wire x1="-17.78" y1="25.4" x2="-17.78" y2="30.48" width="0.1524" layer="91"/>
<wire x1="-17.78" y1="30.48" x2="27.94" y2="30.48" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$23"/>
<pinref part="U$2" gate="G$1" pin="P$28"/>
</segment>
</net>
<net name="N$16" class="0">
<segment>
<wire x1="-25.4" y1="30.48" x2="-20.32" y2="30.48" width="0.1524" layer="91"/>
<wire x1="-20.32" y1="30.48" x2="-20.32" y2="33.02" width="0.1524" layer="91"/>
<wire x1="-20.32" y1="33.02" x2="27.94" y2="33.02" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$22"/>
<pinref part="U$2" gate="G$1" pin="P$30"/>
</segment>
</net>
<net name="N$17" class="0">
<segment>
<wire x1="-25.4" y1="35.56" x2="27.94" y2="35.56" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$21"/>
<pinref part="U$2" gate="G$1" pin="P$32"/>
</segment>
</net>
<net name="N$18" class="0">
<segment>
<wire x1="-53.34" y1="-35.56" x2="-71.12" y2="-35.56" width="0.1524" layer="91"/>
<wire x1="-71.12" y1="-35.56" x2="-71.12" y2="-58.42" width="0.1524" layer="91"/>
<wire x1="-71.12" y1="-58.42" x2="68.58" y2="-58.42" width="0.1524" layer="91"/>
<wire x1="68.58" y1="-58.42" x2="68.58" y2="0" width="0.1524" layer="91"/>
<wire x1="68.58" y1="0" x2="53.34" y2="0" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$2"/>
<pinref part="U$2" gate="G$1" pin="P$3"/>
</segment>
</net>
<net name="N$19" class="0">
<segment>
<wire x1="-53.34" y1="-30.48" x2="-76.2" y2="-30.48" width="0.1524" layer="91"/>
<wire x1="-76.2" y1="-30.48" x2="-76.2" y2="-60.96" width="0.1524" layer="91"/>
<wire x1="-76.2" y1="-60.96" x2="71.12" y2="-60.96" width="0.1524" layer="91"/>
<wire x1="71.12" y1="-60.96" x2="71.12" y2="2.54" width="0.1524" layer="91"/>
<wire x1="71.12" y1="2.54" x2="53.34" y2="2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$3"/>
<pinref part="U$2" gate="G$1" pin="P$5"/>
</segment>
</net>
<net name="N$20" class="0">
<segment>
<wire x1="-53.34" y1="-25.4" x2="-78.74" y2="-25.4" width="0.1524" layer="91"/>
<wire x1="-78.74" y1="-25.4" x2="-78.74" y2="-63.5" width="0.1524" layer="91"/>
<wire x1="-78.74" y1="-63.5" x2="73.66" y2="-63.5" width="0.1524" layer="91"/>
<wire x1="73.66" y1="-63.5" x2="73.66" y2="5.08" width="0.1524" layer="91"/>
<wire x1="73.66" y1="5.08" x2="53.34" y2="5.08" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$4"/>
<pinref part="U$2" gate="G$1" pin="P$7"/>
</segment>
</net>
<net name="N$21" class="0">
<segment>
<wire x1="-53.34" y1="-20.32" x2="-81.28" y2="-20.32" width="0.1524" layer="91"/>
<wire x1="-81.28" y1="-20.32" x2="-81.28" y2="-66.04" width="0.1524" layer="91"/>
<wire x1="-81.28" y1="-66.04" x2="76.2" y2="-66.04" width="0.1524" layer="91"/>
<wire x1="76.2" y1="-66.04" x2="76.2" y2="7.62" width="0.1524" layer="91"/>
<wire x1="76.2" y1="7.62" x2="53.34" y2="7.62" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$5"/>
<pinref part="U$2" gate="G$1" pin="P$9"/>
</segment>
</net>
<net name="N$22" class="0">
<segment>
<wire x1="-53.34" y1="-15.24" x2="-83.82" y2="-15.24" width="0.1524" layer="91"/>
<wire x1="-83.82" y1="-15.24" x2="-83.82" y2="-68.58" width="0.1524" layer="91"/>
<wire x1="-83.82" y1="-68.58" x2="78.74" y2="-68.58" width="0.1524" layer="91"/>
<wire x1="78.74" y1="-68.58" x2="78.74" y2="10.16" width="0.1524" layer="91"/>
<wire x1="78.74" y1="10.16" x2="53.34" y2="10.16" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$6"/>
<pinref part="U$2" gate="G$1" pin="P$11"/>
</segment>
</net>
<net name="N$23" class="0">
<segment>
<wire x1="-53.34" y1="-10.16" x2="-86.36" y2="-10.16" width="0.1524" layer="91"/>
<wire x1="-86.36" y1="-10.16" x2="-86.36" y2="-71.12" width="0.1524" layer="91"/>
<wire x1="-86.36" y1="-71.12" x2="86.36" y2="-71.12" width="0.1524" layer="91"/>
<wire x1="86.36" y1="-71.12" x2="86.36" y2="12.7" width="0.1524" layer="91"/>
<wire x1="86.36" y1="12.7" x2="53.34" y2="12.7" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$7"/>
<pinref part="U$2" gate="G$1" pin="P$13"/>
</segment>
</net>
<net name="N$24" class="0">
<segment>
<wire x1="-53.34" y1="-5.08" x2="-88.9" y2="-5.08" width="0.1524" layer="91"/>
<wire x1="-88.9" y1="-5.08" x2="-88.9" y2="-73.66" width="0.1524" layer="91"/>
<wire x1="-88.9" y1="-73.66" x2="91.44" y2="-73.66" width="0.1524" layer="91"/>
<wire x1="91.44" y1="-73.66" x2="91.44" y2="15.24" width="0.1524" layer="91"/>
<wire x1="91.44" y1="15.24" x2="53.34" y2="15.24" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$8"/>
<pinref part="U$2" gate="G$1" pin="P$15"/>
</segment>
</net>
<net name="N$25" class="0">
<segment>
<wire x1="-53.34" y1="0" x2="-91.44" y2="0" width="0.1524" layer="91"/>
<wire x1="-91.44" y1="0" x2="-91.44" y2="-76.2" width="0.1524" layer="91"/>
<wire x1="-91.44" y1="-76.2" x2="93.98" y2="-76.2" width="0.1524" layer="91"/>
<wire x1="93.98" y1="-76.2" x2="93.98" y2="17.78" width="0.1524" layer="91"/>
<wire x1="93.98" y1="17.78" x2="53.34" y2="17.78" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$9"/>
<pinref part="U$2" gate="G$1" pin="P$17"/>
</segment>
</net>
<net name="N$26" class="0">
<segment>
<wire x1="53.34" y1="20.32" x2="96.52" y2="20.32" width="0.1524" layer="91"/>
<wire x1="96.52" y1="20.32" x2="96.52" y2="-78.74" width="0.1524" layer="91"/>
<wire x1="96.52" y1="-78.74" x2="-93.98" y2="-78.74" width="0.1524" layer="91"/>
<wire x1="-93.98" y1="-78.74" x2="-93.98" y2="5.08" width="0.1524" layer="91"/>
<wire x1="-93.98" y1="5.08" x2="-53.34" y2="5.08" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$10"/>
<pinref part="U$2" gate="G$1" pin="P$19"/>
</segment>
</net>
<net name="N$27" class="0">
<segment>
<wire x1="-53.34" y1="10.16" x2="-96.52" y2="10.16" width="0.1524" layer="91"/>
<wire x1="-96.52" y1="10.16" x2="-96.52" y2="-81.28" width="0.1524" layer="91"/>
<wire x1="-96.52" y1="-81.28" x2="101.6" y2="-81.28" width="0.1524" layer="91"/>
<wire x1="101.6" y1="-81.28" x2="101.6" y2="22.86" width="0.1524" layer="91"/>
<wire x1="101.6" y1="22.86" x2="53.34" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$11"/>
<pinref part="U$2" gate="G$1" pin="P$21"/>
</segment>
</net>
<net name="N$28" class="0">
<segment>
<wire x1="-53.34" y1="15.24" x2="-99.06" y2="15.24" width="0.1524" layer="91"/>
<wire x1="-99.06" y1="15.24" x2="-99.06" y2="-83.82" width="0.1524" layer="91"/>
<wire x1="-99.06" y1="-83.82" x2="106.68" y2="-83.82" width="0.1524" layer="91"/>
<wire x1="106.68" y1="-83.82" x2="106.68" y2="-86.36" width="0.1524" layer="91"/>
<wire x1="106.68" y1="-86.36" x2="109.22" y2="-86.36" width="0.1524" layer="91"/>
<wire x1="109.22" y1="-86.36" x2="109.22" y2="25.4" width="0.1524" layer="91"/>
<wire x1="109.22" y1="25.4" x2="53.34" y2="25.4" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$12"/>
<pinref part="U$2" gate="G$1" pin="P$23"/>
</segment>
</net>
<net name="N$29" class="0">
<segment>
<wire x1="-53.34" y1="20.32" x2="-101.6" y2="20.32" width="0.1524" layer="91"/>
<wire x1="-101.6" y1="20.32" x2="-101.6" y2="-86.36" width="0.1524" layer="91"/>
<wire x1="-101.6" y1="-86.36" x2="101.6" y2="-86.36" width="0.1524" layer="91"/>
<wire x1="101.6" y1="-86.36" x2="101.6" y2="-88.9" width="0.1524" layer="91"/>
<wire x1="101.6" y1="-88.9" x2="111.76" y2="-88.9" width="0.1524" layer="91"/>
<wire x1="111.76" y1="-88.9" x2="111.76" y2="27.94" width="0.1524" layer="91"/>
<wire x1="111.76" y1="27.94" x2="53.34" y2="27.94" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$13"/>
<pinref part="U$2" gate="G$1" pin="P$25"/>
</segment>
</net>
<net name="N$30" class="0">
<segment>
<wire x1="-53.34" y1="25.4" x2="-104.14" y2="25.4" width="0.1524" layer="91"/>
<wire x1="-104.14" y1="25.4" x2="-104.14" y2="-88.9" width="0.1524" layer="91"/>
<wire x1="-104.14" y1="-88.9" x2="99.06" y2="-88.9" width="0.1524" layer="91"/>
<wire x1="99.06" y1="-88.9" x2="99.06" y2="-91.44" width="0.1524" layer="91"/>
<wire x1="99.06" y1="-91.44" x2="114.3" y2="-91.44" width="0.1524" layer="91"/>
<wire x1="114.3" y1="-91.44" x2="114.3" y2="30.48" width="0.1524" layer="91"/>
<wire x1="114.3" y1="30.48" x2="53.34" y2="30.48" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$14"/>
<pinref part="U$2" gate="G$1" pin="P$27"/>
</segment>
</net>
<net name="N$31" class="0">
<segment>
<wire x1="-53.34" y1="30.48" x2="-106.68" y2="30.48" width="0.1524" layer="91"/>
<wire x1="-106.68" y1="30.48" x2="-106.68" y2="-93.98" width="0.1524" layer="91"/>
<wire x1="-106.68" y1="-93.98" x2="116.84" y2="-93.98" width="0.1524" layer="91"/>
<wire x1="116.84" y1="-93.98" x2="116.84" y2="33.02" width="0.1524" layer="91"/>
<wire x1="116.84" y1="33.02" x2="53.34" y2="33.02" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$15"/>
<pinref part="U$2" gate="G$1" pin="P$29"/>
</segment>
</net>
<net name="N$32" class="0">
<segment>
<wire x1="-53.34" y1="35.56" x2="-109.22" y2="35.56" width="0.1524" layer="91"/>
<wire x1="-109.22" y1="35.56" x2="-109.22" y2="-99.06" width="0.1524" layer="91"/>
<wire x1="-109.22" y1="-99.06" x2="119.38" y2="-99.06" width="0.1524" layer="91"/>
<wire x1="119.38" y1="-99.06" x2="119.38" y2="35.56" width="0.1524" layer="91"/>
<wire x1="119.38" y1="35.56" x2="53.34" y2="35.56" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="P$16"/>
<pinref part="U$2" gate="G$1" pin="P$31"/>
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
