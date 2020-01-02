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
<library name="tuckerLib">
<packages>
<package name="SFMC-118-L1-S-D">
<smd name="1" x="-10.795" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="2" x="-10.795" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="3" x="-9.525" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="4" x="-9.525" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="5" x="-8.255" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="6" x="-8.255" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="7" x="-6.985" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="8" x="-6.985" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="9" x="-5.715" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="10" x="-5.715" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="11" x="-4.445" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="12" x="-4.445" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="13" x="-3.175" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="14" x="-3.175" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="15" x="-1.905" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="16" x="-1.905" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="17" x="-0.635" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="18" x="-0.635" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="19" x="0.635" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="20" x="0.635" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="21" x="1.905" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="22" x="1.905" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="23" x="3.175" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="24" x="3.175" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="25" x="4.445" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="26" x="4.445" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="27" x="5.715" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="28" x="5.715" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="29" x="6.985" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="30" x="6.985" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="31" x="8.255" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="32" x="8.255" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="33" x="9.525" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="34" x="9.525" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<smd name="35" x="10.795" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="36" x="10.795" y="-2.032" dx="0.762" dy="2.794" layer="16" rot="R180"/>
<text x="-13.081" y="-2.667" size="1.27" layer="21" ratio="6" rot="SR0">1</text>
<text x="11.938" y="-2.667" size="1.27" layer="21" ratio="6" rot="SR0">35</text>
<wire x1="-11.43" y1="-0.0508" x2="-11.43" y2="-1.7272" width="0.1524" layer="51"/>
<wire x1="-11.43" y1="-1.7272" x2="11.43" y2="-1.7272" width="0.1524" layer="51"/>
<wire x1="11.43" y1="-1.7272" x2="11.43" y2="-0.0508" width="0.1524" layer="51"/>
<wire x1="11.43" y1="-0.0508" x2="-11.43" y2="-0.0508" width="0.1524" layer="51"/>
<text x="-13.081" y="-2.667" size="1.27" layer="51" ratio="6" rot="SR0">1</text>
<text x="11.938" y="-2.667" size="1.27" layer="51" ratio="6" rot="SR0">35</text>
<text x="-11.0744" y="0.3302" size="2.0828" layer="25" ratio="10" rot="SR0">&gt;NAME</text>
<text x="0.4826" y="0.3556" size="2.0828" layer="27" ratio="10" rot="SR0">&gt;VALUE</text>
<text x="-11.43" y="-2.54" size="1.27" layer="22" rot="MR0">2</text>
<text x="13.97" y="-2.54" size="1.27" layer="22" rot="MR0">36</text>
</package>
</packages>
<symbols>
<symbol name="SFMC-118-L1-S-D">
<pin name="1" x="-17.78" y="20.32" length="middle" direction="pas"/>
<pin name="3" x="-17.78" y="17.78" length="middle" direction="pas"/>
<pin name="5" x="-17.78" y="15.24" length="middle" direction="pas"/>
<pin name="7" x="-17.78" y="12.7" length="middle" direction="pas"/>
<pin name="9" x="-17.78" y="10.16" length="middle" direction="pas"/>
<pin name="11" x="-17.78" y="7.62" length="middle" direction="pas"/>
<pin name="13" x="-17.78" y="5.08" length="middle" direction="pas"/>
<pin name="15" x="-17.78" y="2.54" length="middle" direction="pas"/>
<pin name="17" x="-17.78" y="0" length="middle" direction="pas"/>
<pin name="19" x="-17.78" y="-2.54" length="middle" direction="pas"/>
<pin name="21" x="-17.78" y="-5.08" length="middle" direction="pas"/>
<pin name="23" x="-17.78" y="-7.62" length="middle" direction="pas"/>
<pin name="25" x="-17.78" y="-10.16" length="middle" direction="pas"/>
<pin name="27" x="-17.78" y="-12.7" length="middle" direction="pas"/>
<pin name="29" x="-17.78" y="-15.24" length="middle" direction="pas"/>
<pin name="31" x="-17.78" y="-17.78" length="middle" direction="pas"/>
<pin name="33" x="-17.78" y="-20.32" length="middle" direction="pas"/>
<pin name="35" x="-17.78" y="-22.86" length="middle" direction="pas"/>
<pin name="2" x="17.78" y="20.32" length="middle" direction="pas" rot="R180"/>
<pin name="4" x="17.78" y="17.78" length="middle" direction="pas" rot="R180"/>
<pin name="6" x="17.78" y="15.24" length="middle" direction="pas" rot="R180"/>
<pin name="8" x="17.78" y="12.7" length="middle" direction="pas" rot="R180"/>
<pin name="10" x="17.78" y="10.16" length="middle" direction="pas" rot="R180"/>
<pin name="12" x="17.78" y="7.62" length="middle" direction="pas" rot="R180"/>
<pin name="14" x="17.78" y="5.08" length="middle" direction="pas" rot="R180"/>
<pin name="16" x="17.78" y="2.54" length="middle" direction="pas" rot="R180"/>
<pin name="18" x="17.78" y="0" length="middle" direction="pas" rot="R180"/>
<pin name="20" x="17.78" y="-2.54" length="middle" direction="pas" rot="R180"/>
<pin name="22" x="17.78" y="-5.08" length="middle" direction="pas" rot="R180"/>
<pin name="24" x="17.78" y="-7.62" length="middle" direction="pas" rot="R180"/>
<pin name="26" x="17.78" y="-10.16" length="middle" direction="pas" rot="R180"/>
<pin name="28" x="17.78" y="-12.7" length="middle" direction="pas" rot="R180"/>
<pin name="30" x="17.78" y="-15.24" length="middle" direction="pas" rot="R180"/>
<pin name="32" x="17.78" y="-17.78" length="middle" direction="pas" rot="R180"/>
<pin name="34" x="17.78" y="-20.32" length="middle" direction="pas" rot="R180"/>
<pin name="36" x="17.78" y="-22.86" length="middle" direction="pas" rot="R180"/>
<wire x1="-12.7" y1="25.4" x2="-12.7" y2="-27.94" width="0.4064" layer="94"/>
<wire x1="-12.7" y1="-27.94" x2="12.7" y2="-27.94" width="0.4064" layer="94"/>
<wire x1="12.7" y1="-27.94" x2="12.7" y2="25.4" width="0.4064" layer="94"/>
<wire x1="12.7" y1="25.4" x2="-12.7" y2="25.4" width="0.4064" layer="94"/>
<text x="-4.826" y="26.543" size="2.0828" layer="95" ratio="10" rot="SR0">&gt;NAME</text>
<text x="-4.7752" y="-30.8356" size="2.0828" layer="96" ratio="10" rot="SR0">&gt;VALUE</text>
</symbol>
</symbols>
<devicesets>
<deviceset name="SFMC-118-L1-S-D">
<gates>
<gate name="G$1" symbol="SFMC-118-L1-S-D" x="0" y="0"/>
</gates>
<devices>
<device name="" package="SFMC-118-L1-S-D">
<connects>
<connect gate="G$1" pin="1" pad="1"/>
<connect gate="G$1" pin="10" pad="10"/>
<connect gate="G$1" pin="11" pad="11"/>
<connect gate="G$1" pin="12" pad="12"/>
<connect gate="G$1" pin="13" pad="13"/>
<connect gate="G$1" pin="14" pad="14"/>
<connect gate="G$1" pin="15" pad="15"/>
<connect gate="G$1" pin="16" pad="16"/>
<connect gate="G$1" pin="17" pad="17"/>
<connect gate="G$1" pin="18" pad="18"/>
<connect gate="G$1" pin="19" pad="19"/>
<connect gate="G$1" pin="2" pad="2"/>
<connect gate="G$1" pin="20" pad="20"/>
<connect gate="G$1" pin="21" pad="21"/>
<connect gate="G$1" pin="22" pad="22"/>
<connect gate="G$1" pin="23" pad="23"/>
<connect gate="G$1" pin="24" pad="24"/>
<connect gate="G$1" pin="25" pad="25"/>
<connect gate="G$1" pin="26" pad="26"/>
<connect gate="G$1" pin="27" pad="27"/>
<connect gate="G$1" pin="28" pad="28"/>
<connect gate="G$1" pin="29" pad="29"/>
<connect gate="G$1" pin="3" pad="3"/>
<connect gate="G$1" pin="30" pad="30"/>
<connect gate="G$1" pin="31" pad="31"/>
<connect gate="G$1" pin="32" pad="32"/>
<connect gate="G$1" pin="33" pad="33"/>
<connect gate="G$1" pin="34" pad="34"/>
<connect gate="G$1" pin="35" pad="35"/>
<connect gate="G$1" pin="36" pad="36"/>
<connect gate="G$1" pin="4" pad="4"/>
<connect gate="G$1" pin="5" pad="5"/>
<connect gate="G$1" pin="6" pad="6"/>
<connect gate="G$1" pin="7" pad="7"/>
<connect gate="G$1" pin="8" pad="8"/>
<connect gate="G$1" pin="9" pad="9"/>
</connects>
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
</packages>
<symbols>
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
<deviceset name="OMNETICS_NANO_18X2">
<gates>
<gate name="G$1" symbol="OMNETICS_NANO_18X2" x="0" y="0"/>
</gates>
<devices>
<device name="" package="OMNETICS_NANO_18X2">
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
</devices>
</deviceset>
</devicesets>
</library>
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
<part name="U$2" library="Sombeck_lib" deviceset="OMNETICS_NANO_18X2" device=""/>
<part name="GND1" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND2" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND3" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND4" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="AGND1" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="AGND" device=""/>
<part name="AGND2" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="AGND" device=""/>
<part name="U$1" library="tuckerLib" deviceset="SFMC-118-L1-S-D" device=""/>
</parts>
<sheets>
<sheet>
<plain>
</plain>
<instances>
<instance part="U$2" gate="G$1" x="-38.1" y="-2.54"/>
<instance part="GND1" gate="1" x="63.5" y="40.64" rot="R90"/>
<instance part="GND2" gate="1" x="15.24" y="38.1" rot="R270"/>
<instance part="GND3" gate="1" x="-12.7" y="40.64" rot="R90"/>
<instance part="GND4" gate="1" x="-63.5" y="-45.72" rot="R270"/>
<instance part="AGND1" gate="VR1" x="-66.04" y="40.64" rot="R270"/>
<instance part="AGND2" gate="VR1" x="60.96" y="38.1" rot="R90"/>
<instance part="U$1" gate="G$1" x="38.1" y="17.78" rot="R180"/>
</instances>
<busses>
</busses>
<nets>
<net name="GND" class="0">
<segment>
<pinref part="GND1" gate="1" pin="GND"/>
<wire x1="55.88" y1="40.64" x2="60.96" y2="40.64" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="35"/>
</segment>
<segment>
<pinref part="GND2" gate="1" pin="GND"/>
<wire x1="20.32" y1="38.1" x2="17.78" y2="38.1" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="34"/>
</segment>
<segment>
<pinref part="U$2" gate="G$1" pin="GND0"/>
<pinref part="GND4" gate="1" pin="GND"/>
<wire x1="-53.34" y1="-45.72" x2="-60.96" y2="-45.72" width="0.1524" layer="91"/>
</segment>
<segment>
<pinref part="U$2" gate="G$1" pin="GND34"/>
<pinref part="GND3" gate="1" pin="GND"/>
<wire x1="-25.4" y1="40.64" x2="-15.24" y2="40.64" width="0.1524" layer="91"/>
</segment>
</net>
<net name="AGND" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="REF33"/>
<pinref part="AGND1" gate="VR1" pin="AGND"/>
<wire x1="-53.34" y1="40.64" x2="-63.5" y2="40.64" width="0.1524" layer="91"/>
</segment>
<segment>
<pinref part="AGND2" gate="VR1" pin="AGND"/>
<wire x1="55.88" y1="38.1" x2="58.42" y2="38.1" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="33"/>
</segment>
</net>
<net name="N$1" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$1"/>
<wire x1="-53.34" y1="-40.64" x2="-68.58" y2="-40.64" width="0.1524" layer="91"/>
<wire x1="-68.58" y1="-40.64" x2="-68.58" y2="-53.34" width="0.1524" layer="91"/>
<wire x1="-68.58" y1="-53.34" x2="60.96" y2="-53.34" width="0.1524" layer="91"/>
<wire x1="60.96" y1="-53.34" x2="60.96" y2="-2.54" width="0.1524" layer="91"/>
<wire x1="60.96" y1="-2.54" x2="55.88" y2="-2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="1"/>
</segment>
</net>
<net name="N$2" class="0">
<segment>
<wire x1="20.32" y1="-2.54" x2="15.24" y2="-2.54" width="0.1524" layer="91"/>
<wire x1="15.24" y1="-2.54" x2="15.24" y2="-40.64" width="0.1524" layer="91"/>
<pinref part="U$2" gate="G$1" pin="P$2"/>
<wire x1="15.24" y1="-40.64" x2="-25.4" y2="-40.64" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="2"/>
</segment>
</net>
<net name="N$3" class="0">
<segment>
<wire x1="20.32" y1="0" x2="12.7" y2="0" width="0.1524" layer="91"/>
<wire x1="12.7" y1="0" x2="12.7" y2="-35.56" width="0.1524" layer="91"/>
<pinref part="U$2" gate="G$1" pin="P$4"/>
<wire x1="12.7" y1="-35.56" x2="-25.4" y2="-35.56" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="4"/>
</segment>
</net>
<net name="N$4" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$6"/>
<wire x1="-25.4" y1="-30.48" x2="10.16" y2="-30.48" width="0.1524" layer="91"/>
<wire x1="10.16" y1="-30.48" x2="10.16" y2="2.54" width="0.1524" layer="91"/>
<wire x1="10.16" y1="2.54" x2="20.32" y2="2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="6"/>
</segment>
</net>
<net name="N$5" class="0">
<segment>
<wire x1="20.32" y1="5.08" x2="7.62" y2="5.08" width="0.1524" layer="91"/>
<wire x1="7.62" y1="5.08" x2="7.62" y2="-25.4" width="0.1524" layer="91"/>
<pinref part="U$2" gate="G$1" pin="P$8"/>
<wire x1="7.62" y1="-25.4" x2="-25.4" y2="-25.4" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="8"/>
</segment>
</net>
<net name="N$6" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$10"/>
<wire x1="-25.4" y1="-20.32" x2="5.08" y2="-20.32" width="0.1524" layer="91"/>
<wire x1="5.08" y1="-20.32" x2="5.08" y2="7.62" width="0.1524" layer="91"/>
<wire x1="5.08" y1="7.62" x2="20.32" y2="7.62" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="10"/>
</segment>
</net>
<net name="N$7" class="0">
<segment>
<wire x1="20.32" y1="10.16" x2="2.54" y2="10.16" width="0.1524" layer="91"/>
<wire x1="2.54" y1="10.16" x2="2.54" y2="-15.24" width="0.1524" layer="91"/>
<pinref part="U$2" gate="G$1" pin="P$12"/>
<wire x1="2.54" y1="-15.24" x2="-25.4" y2="-15.24" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="12"/>
</segment>
</net>
<net name="N$8" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$14"/>
<wire x1="-25.4" y1="-10.16" x2="0" y2="-10.16" width="0.1524" layer="91"/>
<wire x1="0" y1="-10.16" x2="0" y2="12.7" width="0.1524" layer="91"/>
<wire x1="0" y1="12.7" x2="20.32" y2="12.7" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="14"/>
</segment>
</net>
<net name="N$9" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$16"/>
<wire x1="-25.4" y1="-5.08" x2="-2.54" y2="-5.08" width="0.1524" layer="91"/>
<wire x1="-2.54" y1="-5.08" x2="-2.54" y2="15.24" width="0.1524" layer="91"/>
<wire x1="-2.54" y1="15.24" x2="20.32" y2="15.24" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="16"/>
</segment>
</net>
<net name="N$10" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$18"/>
<wire x1="-25.4" y1="0" x2="-5.08" y2="0" width="0.1524" layer="91"/>
<wire x1="-5.08" y1="0" x2="-5.08" y2="17.78" width="0.1524" layer="91"/>
<wire x1="-5.08" y1="17.78" x2="20.32" y2="17.78" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="18"/>
</segment>
</net>
<net name="N$11" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$20"/>
<wire x1="-25.4" y1="5.08" x2="-7.62" y2="5.08" width="0.1524" layer="91"/>
<wire x1="-7.62" y1="5.08" x2="-7.62" y2="20.32" width="0.1524" layer="91"/>
<wire x1="-7.62" y1="20.32" x2="20.32" y2="20.32" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="20"/>
</segment>
</net>
<net name="N$12" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$22"/>
<wire x1="-25.4" y1="10.16" x2="-10.16" y2="10.16" width="0.1524" layer="91"/>
<wire x1="-10.16" y1="10.16" x2="-10.16" y2="22.86" width="0.1524" layer="91"/>
<wire x1="-10.16" y1="22.86" x2="20.32" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="22"/>
</segment>
</net>
<net name="N$13" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$24"/>
<wire x1="-25.4" y1="15.24" x2="-12.7" y2="15.24" width="0.1524" layer="91"/>
<wire x1="-12.7" y1="15.24" x2="-12.7" y2="25.4" width="0.1524" layer="91"/>
<wire x1="-12.7" y1="25.4" x2="20.32" y2="25.4" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="24"/>
</segment>
</net>
<net name="N$14" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$26"/>
<wire x1="-25.4" y1="20.32" x2="-15.24" y2="20.32" width="0.1524" layer="91"/>
<wire x1="-15.24" y1="20.32" x2="-15.24" y2="27.94" width="0.1524" layer="91"/>
<wire x1="-15.24" y1="27.94" x2="20.32" y2="27.94" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="26"/>
</segment>
</net>
<net name="N$15" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$28"/>
<wire x1="-25.4" y1="25.4" x2="-17.78" y2="25.4" width="0.1524" layer="91"/>
<wire x1="-17.78" y1="25.4" x2="-17.78" y2="30.48" width="0.1524" layer="91"/>
<wire x1="-17.78" y1="30.48" x2="20.32" y2="30.48" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="28"/>
</segment>
</net>
<net name="N$16" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$30"/>
<wire x1="-25.4" y1="30.48" x2="-20.32" y2="30.48" width="0.1524" layer="91"/>
<wire x1="-20.32" y1="30.48" x2="-20.32" y2="33.02" width="0.1524" layer="91"/>
<wire x1="-20.32" y1="33.02" x2="20.32" y2="33.02" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="30"/>
</segment>
</net>
<net name="N$17" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$32"/>
<wire x1="-25.4" y1="35.56" x2="20.32" y2="35.56" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="32"/>
</segment>
</net>
<net name="N$18" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$3"/>
<wire x1="-53.34" y1="-35.56" x2="-71.12" y2="-35.56" width="0.1524" layer="91"/>
<wire x1="-71.12" y1="-35.56" x2="-71.12" y2="-58.42" width="0.1524" layer="91"/>
<wire x1="-71.12" y1="-58.42" x2="68.58" y2="-58.42" width="0.1524" layer="91"/>
<wire x1="68.58" y1="-58.42" x2="68.58" y2="0" width="0.1524" layer="91"/>
<wire x1="68.58" y1="0" x2="55.88" y2="0" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="3"/>
</segment>
</net>
<net name="N$19" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$5"/>
<wire x1="-53.34" y1="-30.48" x2="-76.2" y2="-30.48" width="0.1524" layer="91"/>
<wire x1="-76.2" y1="-30.48" x2="-76.2" y2="-60.96" width="0.1524" layer="91"/>
<wire x1="-76.2" y1="-60.96" x2="71.12" y2="-60.96" width="0.1524" layer="91"/>
<wire x1="71.12" y1="-60.96" x2="71.12" y2="2.54" width="0.1524" layer="91"/>
<wire x1="71.12" y1="2.54" x2="55.88" y2="2.54" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="5"/>
</segment>
</net>
<net name="N$20" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$7"/>
<wire x1="-53.34" y1="-25.4" x2="-78.74" y2="-25.4" width="0.1524" layer="91"/>
<wire x1="-78.74" y1="-25.4" x2="-78.74" y2="-63.5" width="0.1524" layer="91"/>
<wire x1="-78.74" y1="-63.5" x2="73.66" y2="-63.5" width="0.1524" layer="91"/>
<wire x1="73.66" y1="-63.5" x2="73.66" y2="5.08" width="0.1524" layer="91"/>
<wire x1="73.66" y1="5.08" x2="55.88" y2="5.08" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="7"/>
</segment>
</net>
<net name="N$21" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$9"/>
<wire x1="-53.34" y1="-20.32" x2="-81.28" y2="-20.32" width="0.1524" layer="91"/>
<wire x1="-81.28" y1="-20.32" x2="-81.28" y2="-66.04" width="0.1524" layer="91"/>
<wire x1="-81.28" y1="-66.04" x2="76.2" y2="-66.04" width="0.1524" layer="91"/>
<wire x1="76.2" y1="-66.04" x2="76.2" y2="7.62" width="0.1524" layer="91"/>
<wire x1="76.2" y1="7.62" x2="55.88" y2="7.62" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="9"/>
</segment>
</net>
<net name="N$22" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$11"/>
<wire x1="-53.34" y1="-15.24" x2="-83.82" y2="-15.24" width="0.1524" layer="91"/>
<wire x1="-83.82" y1="-15.24" x2="-83.82" y2="-68.58" width="0.1524" layer="91"/>
<wire x1="-83.82" y1="-68.58" x2="78.74" y2="-68.58" width="0.1524" layer="91"/>
<wire x1="78.74" y1="-68.58" x2="78.74" y2="10.16" width="0.1524" layer="91"/>
<wire x1="78.74" y1="10.16" x2="55.88" y2="10.16" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="11"/>
</segment>
</net>
<net name="N$23" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$13"/>
<wire x1="-53.34" y1="-10.16" x2="-86.36" y2="-10.16" width="0.1524" layer="91"/>
<wire x1="-86.36" y1="-10.16" x2="-86.36" y2="-71.12" width="0.1524" layer="91"/>
<wire x1="-86.36" y1="-71.12" x2="86.36" y2="-71.12" width="0.1524" layer="91"/>
<wire x1="86.36" y1="-71.12" x2="86.36" y2="12.7" width="0.1524" layer="91"/>
<wire x1="86.36" y1="12.7" x2="55.88" y2="12.7" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="13"/>
</segment>
</net>
<net name="N$24" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$15"/>
<wire x1="-53.34" y1="-5.08" x2="-88.9" y2="-5.08" width="0.1524" layer="91"/>
<wire x1="-88.9" y1="-5.08" x2="-88.9" y2="-73.66" width="0.1524" layer="91"/>
<wire x1="-88.9" y1="-73.66" x2="91.44" y2="-73.66" width="0.1524" layer="91"/>
<wire x1="91.44" y1="-73.66" x2="91.44" y2="15.24" width="0.1524" layer="91"/>
<wire x1="91.44" y1="15.24" x2="55.88" y2="15.24" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="15"/>
</segment>
</net>
<net name="N$25" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$17"/>
<wire x1="-53.34" y1="0" x2="-91.44" y2="0" width="0.1524" layer="91"/>
<wire x1="-91.44" y1="0" x2="-91.44" y2="-76.2" width="0.1524" layer="91"/>
<wire x1="-91.44" y1="-76.2" x2="93.98" y2="-76.2" width="0.1524" layer="91"/>
<wire x1="93.98" y1="-76.2" x2="93.98" y2="17.78" width="0.1524" layer="91"/>
<wire x1="93.98" y1="17.78" x2="55.88" y2="17.78" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="17"/>
</segment>
</net>
<net name="N$26" class="0">
<segment>
<wire x1="55.88" y1="20.32" x2="96.52" y2="20.32" width="0.1524" layer="91"/>
<wire x1="96.52" y1="20.32" x2="96.52" y2="-78.74" width="0.1524" layer="91"/>
<wire x1="96.52" y1="-78.74" x2="-93.98" y2="-78.74" width="0.1524" layer="91"/>
<wire x1="-93.98" y1="-78.74" x2="-93.98" y2="5.08" width="0.1524" layer="91"/>
<pinref part="U$2" gate="G$1" pin="P$19"/>
<wire x1="-93.98" y1="5.08" x2="-53.34" y2="5.08" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="19"/>
</segment>
</net>
<net name="N$27" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$21"/>
<wire x1="-53.34" y1="10.16" x2="-96.52" y2="10.16" width="0.1524" layer="91"/>
<wire x1="-96.52" y1="10.16" x2="-96.52" y2="-81.28" width="0.1524" layer="91"/>
<wire x1="-96.52" y1="-81.28" x2="101.6" y2="-81.28" width="0.1524" layer="91"/>
<wire x1="101.6" y1="-81.28" x2="101.6" y2="22.86" width="0.1524" layer="91"/>
<wire x1="101.6" y1="22.86" x2="55.88" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="21"/>
</segment>
</net>
<net name="N$28" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$23"/>
<wire x1="-53.34" y1="15.24" x2="-99.06" y2="15.24" width="0.1524" layer="91"/>
<wire x1="-99.06" y1="15.24" x2="-99.06" y2="-83.82" width="0.1524" layer="91"/>
<wire x1="-99.06" y1="-83.82" x2="106.68" y2="-83.82" width="0.1524" layer="91"/>
<wire x1="106.68" y1="-83.82" x2="106.68" y2="-86.36" width="0.1524" layer="91"/>
<wire x1="106.68" y1="-86.36" x2="109.22" y2="-86.36" width="0.1524" layer="91"/>
<wire x1="109.22" y1="-86.36" x2="109.22" y2="25.4" width="0.1524" layer="91"/>
<wire x1="109.22" y1="25.4" x2="55.88" y2="25.4" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="23"/>
</segment>
</net>
<net name="N$29" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$25"/>
<wire x1="-53.34" y1="20.32" x2="-101.6" y2="20.32" width="0.1524" layer="91"/>
<wire x1="-101.6" y1="20.32" x2="-101.6" y2="-86.36" width="0.1524" layer="91"/>
<wire x1="-101.6" y1="-86.36" x2="101.6" y2="-86.36" width="0.1524" layer="91"/>
<wire x1="101.6" y1="-86.36" x2="101.6" y2="-88.9" width="0.1524" layer="91"/>
<wire x1="101.6" y1="-88.9" x2="111.76" y2="-88.9" width="0.1524" layer="91"/>
<wire x1="111.76" y1="-88.9" x2="111.76" y2="27.94" width="0.1524" layer="91"/>
<wire x1="111.76" y1="27.94" x2="55.88" y2="27.94" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="25"/>
</segment>
</net>
<net name="N$30" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$27"/>
<wire x1="-53.34" y1="25.4" x2="-104.14" y2="25.4" width="0.1524" layer="91"/>
<wire x1="-104.14" y1="25.4" x2="-104.14" y2="-88.9" width="0.1524" layer="91"/>
<wire x1="-104.14" y1="-88.9" x2="99.06" y2="-88.9" width="0.1524" layer="91"/>
<wire x1="99.06" y1="-88.9" x2="99.06" y2="-91.44" width="0.1524" layer="91"/>
<wire x1="99.06" y1="-91.44" x2="114.3" y2="-91.44" width="0.1524" layer="91"/>
<wire x1="114.3" y1="-91.44" x2="114.3" y2="30.48" width="0.1524" layer="91"/>
<wire x1="114.3" y1="30.48" x2="55.88" y2="30.48" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="27"/>
</segment>
</net>
<net name="N$31" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$29"/>
<wire x1="-53.34" y1="30.48" x2="-106.68" y2="30.48" width="0.1524" layer="91"/>
<wire x1="-106.68" y1="30.48" x2="-106.68" y2="-93.98" width="0.1524" layer="91"/>
<wire x1="-106.68" y1="-93.98" x2="116.84" y2="-93.98" width="0.1524" layer="91"/>
<wire x1="116.84" y1="-93.98" x2="116.84" y2="33.02" width="0.1524" layer="91"/>
<wire x1="116.84" y1="33.02" x2="55.88" y2="33.02" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="29"/>
</segment>
</net>
<net name="N$32" class="0">
<segment>
<pinref part="U$2" gate="G$1" pin="P$31"/>
<wire x1="-53.34" y1="35.56" x2="-109.22" y2="35.56" width="0.1524" layer="91"/>
<wire x1="-109.22" y1="35.56" x2="-109.22" y2="-99.06" width="0.1524" layer="91"/>
<wire x1="-109.22" y1="-99.06" x2="119.38" y2="-99.06" width="0.1524" layer="91"/>
<wire x1="119.38" y1="-99.06" x2="119.38" y2="35.56" width="0.1524" layer="91"/>
<wire x1="119.38" y1="35.56" x2="55.88" y2="35.56" width="0.1524" layer="91"/>
<pinref part="U$1" gate="G$1" pin="31"/>
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
