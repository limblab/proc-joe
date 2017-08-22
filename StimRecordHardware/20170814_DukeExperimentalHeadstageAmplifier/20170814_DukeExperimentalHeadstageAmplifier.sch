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
<library name="tuckerLib2">
<packages>
<package name="FTSH-118-01-L-MT-INVERTED">
<smd name="1" x="-10.795" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="3" x="-9.525" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="5" x="-8.255" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="7" x="-6.985" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="9" x="-5.715" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="11" x="-4.445" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="13" x="-3.175" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="15" x="-1.905" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<wire x1="11.43" y1="-1.7018" x2="11.43" y2="-0.8382" width="0.1524" layer="21"/>
<wire x1="11.43" y1="-2.9972" x2="11.43" y2="0.4572" width="0.1524" layer="51"/>
<text x="-5.588" y="-5.2832" size="2.0828" layer="25" ratio="10" rot="SR0">&gt;NAME</text>
<pad name="2" x="-10.795" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="4" x="-9.525" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="6" x="-8.255" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="8" x="-6.985" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="10" x="-5.715" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="12" x="-4.445" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="14" x="-3.175" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="16" x="-1.905" y="-2.54" drill="0.6604" rot="R180"/>
<smd name="17" x="-0.635" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="19" x="0.635" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="21" x="1.905" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="23" x="3.175" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="25" x="4.445" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="27" x="5.715" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="29" x="6.985" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="31" x="8.255" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<wire x1="11.43" y1="0.4572" x2="-11.43" y2="0.4572" width="0.1524" layer="51"/>
<wire x1="-11.43" y1="-2.9972" x2="11.43" y2="-2.9972" width="0.1524" layer="51"/>
<pad name="18" x="-0.635" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="20" x="0.635" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="22" x="1.905" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="24" x="3.175" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="26" x="4.445" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="28" x="5.715" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="30" x="6.985" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="32" x="8.255" y="-2.54" drill="0.6604" rot="R180"/>
<smd name="33" x="9.525" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<smd name="35" x="10.795" y="0.762" dx="0.762" dy="2.794" layer="1" rot="R180"/>
<wire x1="-11.43" y1="-0.8382" x2="-11.43" y2="-1.7018" width="0.1524" layer="21"/>
<wire x1="-11.43" y1="0.4572" x2="-11.43" y2="-2.9972" width="0.1524" layer="51"/>
<pad name="34" x="9.525" y="-2.54" drill="0.6604" rot="R180"/>
<pad name="36" x="10.795" y="-2.54" drill="0.6604" rot="R180"/>
<text x="-13.462" y="0.127" size="1.27" layer="21">01</text>
<text x="-13.589" y="-2.794" size="1.27" layer="21">02</text>
</package>
<package name="FTSH-118-01-L-MT">
<smd name="1" x="-9.525" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="3" x="-8.255" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="5" x="-6.985" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="7" x="-5.715" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="9" x="-4.445" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="11" x="-3.175" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="13" x="-1.905" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="15" x="-0.635" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<wire x1="-10.16" y1="0.4318" x2="-10.16" y2="-0.4318" width="0.1524" layer="21"/>
<text x="-11.176" y="0.127" size="1.27" layer="21" ratio="6" rot="SR0">2</text>
<text x="-11.176" y="-2.667" size="1.27" layer="21" ratio="6" rot="SR0">1</text>
<wire x1="-10.16" y1="1.7272" x2="-10.16" y2="-1.7272" width="0.1524" layer="51"/>
<text x="-11.176" y="0.127" size="1.27" layer="51" ratio="6" rot="SR0">2</text>
<text x="-11.176" y="-2.667" size="1.27" layer="51" ratio="6" rot="SR0">1</text>
<text x="-10.287" y="-5.6642" size="2.0828" layer="25" ratio="10" rot="SR0">&gt;NAME</text>
<pad name="2" x="-9.525" y="1.27" drill="0.6604"/>
<pad name="4" x="-8.255" y="1.27" drill="0.6604"/>
<pad name="6" x="-6.985" y="1.27" drill="0.6604"/>
<pad name="8" x="-5.715" y="1.27" drill="0.6604"/>
<pad name="10" x="-4.445" y="1.27" drill="0.6604"/>
<pad name="12" x="-3.175" y="1.27" drill="0.6604"/>
<pad name="14" x="-1.905" y="1.27" drill="0.6604"/>
<pad name="16" x="-0.635" y="1.27" drill="0.6604"/>
<smd name="17" x="0.635" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="19" x="1.905" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="21" x="3.175" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="23" x="4.445" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="25" x="5.715" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="27" x="6.985" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="29" x="8.255" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="31" x="9.525" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<wire x1="-10.16" y1="-1.7272" x2="12.7" y2="-1.7272" width="0.1524" layer="51"/>
<wire x1="12.7" y1="1.7272" x2="-10.16" y2="1.7272" width="0.1524" layer="51"/>
<pad name="18" x="0.635" y="1.27" drill="0.6604"/>
<pad name="20" x="1.905" y="1.27" drill="0.6604"/>
<pad name="22" x="3.175" y="1.27" drill="0.6604"/>
<pad name="24" x="4.445" y="1.27" drill="0.6604"/>
<pad name="26" x="5.715" y="1.27" drill="0.6604"/>
<pad name="28" x="6.985" y="1.27" drill="0.6604"/>
<pad name="30" x="8.255" y="1.27" drill="0.6604"/>
<pad name="32" x="9.525" y="1.27" drill="0.6604"/>
<smd name="33" x="10.795" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="35" x="12.065" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<wire x1="12.7" y1="-0.4318" x2="12.7" y2="0.4318" width="0.1524" layer="21"/>
<text x="12.573" y="0.127" size="1.27" layer="21" ratio="6" rot="SR0">36</text>
<text x="12.573" y="-2.667" size="1.27" layer="21" ratio="6" rot="SR0">35</text>
<wire x1="12.7" y1="-1.7272" x2="12.7" y2="1.7272" width="0.1524" layer="51"/>
<text x="12.573" y="0.127" size="1.27" layer="51" ratio="6" rot="SR0">36</text>
<text x="12.573" y="-2.667" size="1.27" layer="51" ratio="6" rot="SR0">35</text>
<pad name="34" x="10.795" y="1.27" drill="0.6604"/>
<pad name="36" x="12.065" y="1.27" drill="0.6604"/>
</package>
<package name="FTSH-132-01-DV">
<smd name="1" x="-10.795" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="2" x="-10.795" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="3" x="-9.525" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="4" x="-9.525" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="5" x="-8.255" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="6" x="-8.255" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="7" x="-6.985" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="8" x="-6.985" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="9" x="-5.715" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="10" x="-5.715" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="11" x="-4.445" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="12" x="-4.445" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="13" x="-3.175" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="14" x="-3.175" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="15" x="-1.905" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="16" x="-1.905" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="17" x="-0.635" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="18" x="-0.635" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="19" x="0.635" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="20" x="0.635" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="21" x="1.905" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="22" x="1.905" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="23" x="3.175" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="24" x="3.175" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="25" x="4.445" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="26" x="4.445" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="27" x="5.715" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="28" x="5.715" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="29" x="6.985" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="30" x="6.985" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="31" x="8.255" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="32" x="8.255" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="33" x="9.525" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="34" x="9.525" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="35" x="10.795" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="36" x="10.795" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<wire x1="-11.43" y1="0.4318" x2="-11.43" y2="-0.4318" width="0.1524" layer="21"/>
<wire x1="-11.43" y1="1.7272" x2="-11.43" y2="-1.7272" width="0.1524" layer="21"/>
<wire x1="29.1846" y1="1.778" x2="-11.43" y2="1.7272" width="0.1524" layer="21"/>
<text x="-5.2324" y="4.9022" size="2.0828" layer="25" ratio="10" rot="SR0">&gt;NAME</text>
<text x="-5.8674" y="-6.7564" size="2.0828" layer="27" ratio="10" rot="SR0">&gt;VALUE</text>
<smd name="37" x="12.065" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="39" x="13.335" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="41" x="14.605" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="43" x="15.875" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="45" x="17.145" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="47" x="18.415" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="49" x="19.685" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="51" x="20.955" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="53" x="22.225" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="55" x="23.495" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="57" x="24.765" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="59" x="26.035" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="61" x="27.305" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="63" x="28.575" y="-2.032" dx="0.762" dy="2.794" layer="1"/>
<wire x1="29.21" y1="1.8034" x2="29.21" y2="-1.778" width="0.127" layer="21"/>
<smd name="38" x="12.065" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="40" x="13.335" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="42" x="14.605" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="44" x="15.875" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="46" x="17.145" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="48" x="18.415" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="50" x="19.685" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="52" x="20.955" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="54" x="22.225" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="56" x="23.495" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="58" x="24.765" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="60" x="26.035" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="62" x="27.305" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<smd name="64" x="28.575" y="2.032" dx="0.762" dy="2.794" layer="1"/>
<wire x1="-11.43" y1="-1.778" x2="29.21" y2="-1.778" width="0.127" layer="21"/>
</package>
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
<symbol name="FTSH-118-01-L-MT">
<pin name="1" x="-17.78" y="22.86" length="middle" direction="pas"/>
<pin name="3" x="-17.78" y="20.32" length="middle" direction="pas"/>
<pin name="5" x="-17.78" y="17.78" length="middle" direction="pas"/>
<pin name="7" x="-17.78" y="15.24" length="middle" direction="pas"/>
<pin name="9" x="-17.78" y="12.7" length="middle" direction="pas"/>
<pin name="11" x="-17.78" y="10.16" length="middle" direction="pas"/>
<pin name="13" x="-17.78" y="7.62" length="middle" direction="pas"/>
<pin name="15" x="-17.78" y="5.08" length="middle" direction="pas"/>
<pin name="2" x="17.78" y="22.86" length="middle" direction="pas" rot="R180"/>
<pin name="4" x="17.78" y="20.32" length="middle" direction="pas" rot="R180"/>
<pin name="6" x="17.78" y="17.78" length="middle" direction="pas" rot="R180"/>
<pin name="8" x="17.78" y="15.24" length="middle" direction="pas" rot="R180"/>
<pin name="10" x="17.78" y="12.7" length="middle" direction="pas" rot="R180"/>
<pin name="12" x="17.78" y="10.16" length="middle" direction="pas" rot="R180"/>
<pin name="14" x="17.78" y="7.62" length="middle" direction="pas" rot="R180"/>
<pin name="16" x="17.78" y="5.08" length="middle" direction="pas" rot="R180"/>
<wire x1="-12.7" y1="27.94" x2="-12.7" y2="-25.4" width="0.4064" layer="94"/>
<wire x1="12.7" y1="-25.4" x2="12.7" y2="27.94" width="0.4064" layer="94"/>
<wire x1="12.7" y1="27.94" x2="-12.7" y2="27.94" width="0.4064" layer="94"/>
<text x="-5.0038" y="29.8958" size="2.0828" layer="95" ratio="10" rot="SR0">&gt;NAME</text>
<pin name="17" x="-17.78" y="2.54" length="middle" direction="pas"/>
<pin name="18" x="17.78" y="2.54" length="middle" direction="pas" rot="R180"/>
<pin name="19" x="-17.78" y="0" length="middle" direction="pas"/>
<pin name="20" x="17.78" y="0" length="middle" direction="pas" rot="R180"/>
<pin name="21" x="-17.78" y="-2.54" length="middle" direction="pas"/>
<pin name="22" x="17.78" y="-2.54" length="middle" direction="pas" rot="R180"/>
<pin name="23" x="-17.78" y="-5.08" length="middle" direction="pas"/>
<pin name="24" x="17.78" y="-5.08" length="middle" direction="pas" rot="R180"/>
<pin name="25" x="-17.78" y="-7.62" length="middle" direction="pas"/>
<pin name="26" x="17.78" y="-7.62" length="middle" direction="pas" rot="R180"/>
<pin name="27" x="-17.78" y="-10.16" length="middle" direction="pas"/>
<pin name="28" x="17.78" y="-10.16" length="middle" direction="pas" rot="R180"/>
<pin name="29" x="-17.78" y="-12.7" length="middle" direction="pas"/>
<pin name="30" x="17.78" y="-12.7" length="middle" direction="pas" rot="R180"/>
<pin name="31" x="-17.78" y="-15.24" length="middle" direction="pas"/>
<pin name="32" x="17.78" y="-15.24" length="middle" direction="pas" rot="R180"/>
<wire x1="-12.7" y1="-25.4" x2="12.7" y2="-25.4" width="0.4064" layer="94"/>
<text x="-6.8834" y="-28.1178" size="2.0828" layer="96" ratio="10" rot="SR0">&gt;VALUE</text>
<pin name="33" x="-17.78" y="-17.78" length="middle" direction="pas"/>
<pin name="34" x="17.78" y="-17.78" length="middle" direction="pas" rot="R180"/>
<pin name="35" x="-17.78" y="-20.32" length="middle" direction="pas"/>
<pin name="36" x="17.78" y="-20.32" length="middle" direction="pas" rot="R180"/>
</symbol>
<symbol name="FTSH-132-01-L-DV">
<pin name="1" x="-17.78" y="30.48" length="middle" direction="pas"/>
<pin name="3" x="-17.78" y="27.94" length="middle" direction="pas"/>
<pin name="5" x="-17.78" y="25.4" length="middle" direction="pas"/>
<pin name="7" x="-17.78" y="22.86" length="middle" direction="pas"/>
<pin name="9" x="-17.78" y="20.32" length="middle" direction="pas"/>
<pin name="11" x="-17.78" y="17.78" length="middle" direction="pas"/>
<pin name="13" x="-17.78" y="15.24" length="middle" direction="pas"/>
<pin name="15" x="-17.78" y="12.7" length="middle" direction="pas"/>
<pin name="17" x="-17.78" y="10.16" length="middle" direction="pas"/>
<pin name="19" x="-17.78" y="7.62" length="middle" direction="pas"/>
<pin name="21" x="-17.78" y="5.08" length="middle" direction="pas"/>
<pin name="23" x="-17.78" y="2.54" length="middle" direction="pas"/>
<pin name="25" x="-17.78" y="0" length="middle" direction="pas"/>
<pin name="27" x="-17.78" y="-2.54" length="middle" direction="pas"/>
<pin name="29" x="-17.78" y="-5.08" length="middle" direction="pas"/>
<pin name="31" x="-17.78" y="-7.62" length="middle" direction="pas"/>
<pin name="33" x="-17.78" y="-10.16" length="middle" direction="pas"/>
<pin name="35" x="-17.78" y="-12.7" length="middle" direction="pas"/>
<pin name="2" x="17.78" y="30.48" length="middle" direction="pas" rot="R180"/>
<pin name="4" x="17.78" y="27.94" length="middle" direction="pas" rot="R180"/>
<pin name="6" x="17.78" y="25.4" length="middle" direction="pas" rot="R180"/>
<pin name="8" x="17.78" y="22.86" length="middle" direction="pas" rot="R180"/>
<pin name="10" x="17.78" y="20.32" length="middle" direction="pas" rot="R180"/>
<pin name="12" x="17.78" y="17.78" length="middle" direction="pas" rot="R180"/>
<pin name="14" x="17.78" y="15.24" length="middle" direction="pas" rot="R180"/>
<pin name="16" x="17.78" y="12.7" length="middle" direction="pas" rot="R180"/>
<pin name="18" x="17.78" y="10.16" length="middle" direction="pas" rot="R180"/>
<pin name="20" x="17.78" y="7.62" length="middle" direction="pas" rot="R180"/>
<pin name="22" x="17.78" y="5.08" length="middle" direction="pas" rot="R180"/>
<pin name="24" x="17.78" y="2.54" length="middle" direction="pas" rot="R180"/>
<pin name="26" x="17.78" y="0" length="middle" direction="pas" rot="R180"/>
<pin name="28" x="17.78" y="-2.54" length="middle" direction="pas" rot="R180"/>
<pin name="30" x="17.78" y="-5.08" length="middle" direction="pas" rot="R180"/>
<pin name="32" x="17.78" y="-7.62" length="middle" direction="pas" rot="R180"/>
<pin name="34" x="17.78" y="-10.16" length="middle" direction="pas" rot="R180"/>
<pin name="36" x="17.78" y="-12.7" length="middle" direction="pas" rot="R180"/>
<wire x1="-12.7" y1="35.56" x2="-12.7" y2="-50.8" width="0.4064" layer="94"/>
<wire x1="-12.7" y1="-50.8" x2="12.7" y2="-50.8" width="0.4064" layer="94"/>
<wire x1="12.7" y1="-50.8" x2="12.7" y2="35.56" width="0.4064" layer="94"/>
<wire x1="12.7" y1="35.56" x2="-12.7" y2="35.56" width="0.4064" layer="94"/>
<text x="-4.826" y="36.703" size="2.0828" layer="95" ratio="10" rot="SR0">&gt;NAME</text>
<text x="-7.3152" y="-53.6956" size="2.0828" layer="96" ratio="10" rot="SR0">&gt;VALUE</text>
<pin name="37" x="-17.78" y="-15.24" length="middle"/>
<pin name="39" x="-17.78" y="-17.78" length="middle"/>
<pin name="41" x="-17.78" y="-20.32" length="middle"/>
<pin name="43" x="-17.78" y="-22.86" length="middle"/>
<pin name="45" x="-17.78" y="-25.4" length="middle"/>
<pin name="47" x="-17.78" y="-27.94" length="middle"/>
<pin name="49" x="-17.78" y="-30.48" length="middle"/>
<pin name="51" x="-17.78" y="-33.02" length="middle"/>
<pin name="53" x="-17.78" y="-35.56" length="middle"/>
<pin name="55" x="-17.78" y="-38.1" length="middle"/>
<pin name="57" x="-17.78" y="-40.64" length="middle"/>
<pin name="59" x="-17.78" y="-43.18" length="middle"/>
<pin name="61" x="-17.78" y="-45.72" length="middle"/>
<pin name="63" x="-17.78" y="-48.26" length="middle"/>
<pin name="38" x="17.78" y="-15.24" length="middle" rot="R180"/>
<pin name="40" x="17.78" y="-17.78" length="middle" rot="R180"/>
<pin name="42" x="17.78" y="-20.32" length="middle" rot="R180"/>
<pin name="44" x="17.78" y="-22.86" length="middle" rot="R180"/>
<pin name="46" x="17.78" y="-25.4" length="middle" rot="R180"/>
<pin name="48" x="17.78" y="-27.94" length="middle" rot="R180"/>
<pin name="50" x="17.78" y="-30.48" length="middle" rot="R180"/>
<pin name="52" x="17.78" y="-33.02" length="middle" rot="R180"/>
<pin name="54" x="17.78" y="-35.56" length="middle" rot="R180"/>
<pin name="56" x="17.78" y="-38.1" length="middle" rot="R180"/>
<pin name="58" x="17.78" y="-40.64" length="middle" rot="R180"/>
<pin name="60" x="17.78" y="-43.18" length="middle" rot="R180"/>
<pin name="62" x="17.78" y="-45.72" length="middle" rot="R180"/>
<pin name="64" x="17.78" y="-48.26" length="middle" rot="R180"/>
</symbol>
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
<deviceset name="FTSH-118-01-L-MT">
<gates>
<gate name="G$1" symbol="FTSH-118-01-L-MT" x="0" y="0"/>
</gates>
<devices>
<device name="REGULAR" package="FTSH-118-01-L-MT">
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
<device name="INVERTED" package="FTSH-118-01-L-MT-INVERTED">
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
<deviceset name="FTSH-132-01-L-DV">
<gates>
<gate name="G$1" symbol="FTSH-132-01-L-DV" x="0" y="7.62"/>
</gates>
<devices>
<device name="" package="FTSH-132-01-DV">
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
<connect gate="G$1" pin="37" pad="37"/>
<connect gate="G$1" pin="38" pad="38"/>
<connect gate="G$1" pin="39" pad="39"/>
<connect gate="G$1" pin="4" pad="4"/>
<connect gate="G$1" pin="40" pad="40"/>
<connect gate="G$1" pin="41" pad="41"/>
<connect gate="G$1" pin="42" pad="42"/>
<connect gate="G$1" pin="43" pad="43"/>
<connect gate="G$1" pin="44" pad="44"/>
<connect gate="G$1" pin="45" pad="45"/>
<connect gate="G$1" pin="46" pad="46"/>
<connect gate="G$1" pin="47" pad="47"/>
<connect gate="G$1" pin="48" pad="48"/>
<connect gate="G$1" pin="49" pad="49"/>
<connect gate="G$1" pin="5" pad="5"/>
<connect gate="G$1" pin="50" pad="50"/>
<connect gate="G$1" pin="51" pad="51"/>
<connect gate="G$1" pin="52" pad="52"/>
<connect gate="G$1" pin="53" pad="53"/>
<connect gate="G$1" pin="54" pad="54"/>
<connect gate="G$1" pin="55" pad="55"/>
<connect gate="G$1" pin="56" pad="56"/>
<connect gate="G$1" pin="57" pad="57"/>
<connect gate="G$1" pin="58" pad="58"/>
<connect gate="G$1" pin="59" pad="59"/>
<connect gate="G$1" pin="6" pad="6"/>
<connect gate="G$1" pin="60" pad="60"/>
<connect gate="G$1" pin="61" pad="61"/>
<connect gate="G$1" pin="62" pad="62"/>
<connect gate="G$1" pin="63" pad="63"/>
<connect gate="G$1" pin="64" pad="64"/>
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
<symbol name="0V" library_version="1">
<wire x1="-1.905" y1="0" x2="1.905" y2="0" width="0.254" layer="94"/>
<text x="-1.905" y="-2.54" size="1.778" layer="96">&gt;VALUE</text>
<pin name="0V" x="0" y="2.54" visible="off" length="short" direction="sup" rot="R270"/>
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
<deviceset name="0V" prefix="GND" library_version="1">
<description>&lt;b&gt;SUPPLY SYMBOL&lt;/b&gt;</description>
<gates>
<gate name="1" symbol="0V" x="0" y="0"/>
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
<library name="con-molex" urn="urn:adsk.eagle:library:165">
<description>&lt;b&gt;Molex Connectors&lt;/b&gt;&lt;p&gt;
&lt;author&gt;Created by librarian@cadsoft.de&lt;/author&gt;</description>
<packages>
<package name="22-23-2021" library_version="1">
<description>.100" (2.54mm) Center Headers - 2 Pin</description>
<wire x1="-2.54" y1="3.175" x2="2.54" y2="3.175" width="0.254" layer="21"/>
<wire x1="2.54" y1="3.175" x2="2.54" y2="1.27" width="0.254" layer="21"/>
<wire x1="2.54" y1="1.27" x2="2.54" y2="-3.175" width="0.254" layer="21"/>
<wire x1="2.54" y1="-3.175" x2="-2.54" y2="-3.175" width="0.254" layer="21"/>
<wire x1="-2.54" y1="-3.175" x2="-2.54" y2="1.27" width="0.254" layer="21"/>
<wire x1="-2.54" y1="1.27" x2="-2.54" y2="3.175" width="0.254" layer="21"/>
<wire x1="-2.54" y1="1.27" x2="2.54" y2="1.27" width="0.254" layer="21"/>
<pad name="1" x="-1.27" y="0" drill="1" shape="long" rot="R90"/>
<pad name="2" x="1.27" y="0" drill="1" shape="long" rot="R90"/>
<text x="-2.54" y="3.81" size="1.016" layer="25" ratio="10">&gt;NAME</text>
<text x="-2.54" y="-5.08" size="1.016" layer="27" ratio="10">&gt;VALUE</text>
</package>
</packages>
<symbols>
<symbol name="MV" library_version="1">
<wire x1="1.27" y1="0" x2="0" y2="0" width="0.6096" layer="94"/>
<text x="2.54" y="-0.762" size="1.524" layer="95">&gt;NAME</text>
<text x="-0.762" y="1.397" size="1.778" layer="96">&gt;VALUE</text>
<pin name="S" x="-2.54" y="0" visible="off" length="short" direction="pas"/>
</symbol>
<symbol name="M" library_version="1">
<wire x1="1.27" y1="0" x2="0" y2="0" width="0.6096" layer="94"/>
<text x="2.54" y="-0.762" size="1.524" layer="95">&gt;NAME</text>
<pin name="S" x="-2.54" y="0" visible="off" length="short" direction="pas"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="22-23-2021" prefix="X" library_version="1">
<description>.100" (2.54mm) Center Header - 2 Pin</description>
<gates>
<gate name="-1" symbol="MV" x="0" y="0" addlevel="always" swaplevel="1"/>
<gate name="-2" symbol="M" x="0" y="-2.54" addlevel="always" swaplevel="1"/>
</gates>
<devices>
<device name="" package="22-23-2021">
<connects>
<connect gate="-1" pin="S" pad="1"/>
<connect gate="-2" pin="S" pad="2"/>
</connects>
<technologies>
<technology name="">
<attribute name="MF" value="MOLEX" constant="no"/>
<attribute name="MPN" value="22-23-2021" constant="no"/>
<attribute name="OC_FARNELL" value="1462926" constant="no"/>
<attribute name="OC_NEWARK" value="25C3832" constant="no"/>
</technology>
</technologies>
</device>
</devices>
</deviceset>
</devicesets>
</library>
<library name="Sombeck_lib">
<packages>
<package name="VIA">
<pad name="P$1" x="0" y="0" drill="0.5" diameter="0.8128"/>
</package>
<package name="LPPB051NFFN-RC">
<pad name="P$1" x="0" y="0" drill="0.5"/>
<pad name="P$2" x="-1.27" y="0" drill="0.5"/>
<pad name="P$3" x="-2.54" y="0" drill="0.5"/>
<pad name="P$4" x="-3.81" y="0" drill="0.5"/>
<pad name="P$5" x="-5.08" y="0" drill="0.5"/>
</package>
</packages>
<symbols>
<symbol name="VIA">
<pin name="P$1" x="-5.08" y="0" length="middle"/>
</symbol>
<symbol name="LPPB051NFFN-RC">
<pin name="P$1" x="-5.08" y="5.08" length="middle"/>
<pin name="P$2" x="-5.08" y="2.54" length="middle"/>
<pin name="P$3" x="-5.08" y="0" length="middle"/>
<pin name="P$4" x="-5.08" y="-2.54" length="middle"/>
<pin name="P$5" x="-5.08" y="-5.08" length="middle"/>
<wire x1="-2.54" y1="7.62" x2="-2.54" y2="-7.62" width="0.254" layer="94"/>
<wire x1="-2.54" y1="-7.62" x2="2.54" y2="-7.62" width="0.254" layer="94"/>
<wire x1="2.54" y1="-7.62" x2="2.54" y2="7.62" width="0.254" layer="94"/>
<wire x1="2.54" y1="7.62" x2="-2.54" y2="7.62" width="0.254" layer="94"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="VIA">
<gates>
<gate name="G$1" symbol="VIA" x="0" y="0"/>
</gates>
<devices>
<device name="" package="VIA">
<connects>
<connect gate="G$1" pin="P$1" pad="P$1"/>
</connects>
<technologies>
<technology name=""/>
</technologies>
</device>
</devices>
</deviceset>
<deviceset name="LPPB051NFFN-RC">
<gates>
<gate name="G$1" symbol="LPPB051NFFN-RC" x="0" y="0"/>
</gates>
<devices>
<device name="" package="LPPB051NFFN-RC">
<connects>
<connect gate="G$1" pin="P$1" pad="P$1"/>
<connect gate="G$1" pin="P$2" pad="P$2"/>
<connect gate="G$1" pin="P$3" pad="P$3"/>
<connect gate="G$1" pin="P$4" pad="P$4"/>
<connect gate="G$1" pin="P$5" pad="P$5"/>
</connects>
<technologies>
<technology name=""/>
</technologies>
</device>
</devices>
</deviceset>
</devicesets>
</library>
<library name="supply2" urn="urn:adsk.eagle:library:372">
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
<symbol name="0V" library_version="1">
<wire x1="-1.27" y1="0" x2="1.27" y2="0" width="0.4064" layer="94"/>
<text x="-1.143" y="-2.032" size="1.778" layer="96">&gt;VALUE</text>
<pin name="0V" x="0" y="2.54" visible="off" length="short" direction="sup" rot="R270"/>
</symbol>
</symbols>
<devicesets>
<deviceset name="0V" prefix="SUPPLY" library_version="1">
<description>&lt;b&gt;SUPPLY SYMBOL&lt;/b&gt;</description>
<gates>
<gate name="0V" symbol="0V" x="0" y="0"/>
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
<part name="GND1" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND2" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="DUKE-STIM" library="con-molex" library_urn="urn:adsk.eagle:library:165" deviceset="22-23-2021" device=""/>
<part name="DUKE-REC" library="con-molex" library_urn="urn:adsk.eagle:library:165" deviceset="22-23-2021" device=""/>
<part name="GND3" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND4" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND7" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="0V" device=""/>
<part name="STIM" library="tuckerLib2" deviceset="FTSH-118-01-L-MT" device="REGULAR"/>
<part name="BUFFER" library="tuckerLib2" deviceset="FTSH-118-01-L-MT" device="REGULAR"/>
<part name="U$3" library="tuckerLib2" deviceset="FTSH-132-01-L-DV" device=""/>
<part name="GND5" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="GND6" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="U$1" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$2" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$4" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$5" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$6" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$7" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$8" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$9" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$10" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$11" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$12" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$13" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$14" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$15" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$16" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$17" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$18" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$19" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$20" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$21" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$22" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$23" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$24" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$25" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$26" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$27" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$28" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$29" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$30" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$31" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$32" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="U$33" library="Sombeck_lib" deviceset="VIA" device=""/>
<part name="GND9" library="supply1" library_urn="urn:adsk.eagle:library:371" deviceset="GND" device=""/>
<part name="SUPPLY1" library="supply2" library_urn="urn:adsk.eagle:library:372" deviceset="0V" device=""/>
<part name="U$34" library="Sombeck_lib" deviceset="LPPB051NFFN-RC" device=""/>
<part name="HEADSTAGE" library="tuckerLib2" deviceset="SFMC-118-L1-S-D" device=""/>
<part name="U$35" library="Sombeck_lib" deviceset="LPPB051NFFN-RC" device=""/>
<part name="U$36" library="Sombeck_lib" deviceset="LPPB051NFFN-RC" device=""/>
<part name="SUPPLY2" library="supply2" library_urn="urn:adsk.eagle:library:372" deviceset="0V" device=""/>
</parts>
<sheets>
<sheet>
<plain>
</plain>
<instances>
<instance part="GND1" gate="1" x="-127" y="-35.56" rot="R270"/>
<instance part="GND2" gate="1" x="-78.74" y="-33.02" rot="R90"/>
<instance part="DUKE-STIM" gate="-1" x="104.14" y="58.42"/>
<instance part="DUKE-STIM" gate="-2" x="104.14" y="55.88"/>
<instance part="DUKE-REC" gate="-1" x="101.6" y="76.2"/>
<instance part="DUKE-REC" gate="-2" x="101.6" y="73.66"/>
<instance part="GND3" gate="1" x="35.56" y="-17.78" rot="R90"/>
<instance part="GND4" gate="1" x="35.56" y="-58.42" rot="R90"/>
<instance part="GND7" gate="1" x="35.56" y="-10.16" rot="R90"/>
<instance part="STIM" gate="G$1" x="-104.14" y="-15.24"/>
<instance part="BUFFER" gate="G$1" x="-15.24" y="129.54" rot="R90"/>
<instance part="U$3" gate="G$1" x="-7.62" y="40.64" rot="R90"/>
<instance part="GND5" gate="1" x="5.08" y="157.48" rot="R180"/>
<instance part="GND6" gate="1" x="27.94" y="111.76" rot="R90"/>
<instance part="U$1" gate="G$1" x="-78.74" y="7.62"/>
<instance part="U$2" gate="G$1" x="-78.74" y="5.08"/>
<instance part="U$4" gate="G$1" x="-78.74" y="2.54"/>
<instance part="U$5" gate="G$1" x="-78.74" y="0"/>
<instance part="U$6" gate="G$1" x="-78.74" y="-2.54"/>
<instance part="U$7" gate="G$1" x="-78.74" y="-5.08"/>
<instance part="U$8" gate="G$1" x="-78.74" y="-7.62"/>
<instance part="U$9" gate="G$1" x="-78.74" y="-10.16"/>
<instance part="U$10" gate="G$1" x="-78.74" y="-12.7"/>
<instance part="U$11" gate="G$1" x="-78.74" y="-15.24"/>
<instance part="U$12" gate="G$1" x="-78.74" y="-17.78"/>
<instance part="U$13" gate="G$1" x="-78.74" y="-20.32"/>
<instance part="U$14" gate="G$1" x="-78.74" y="-22.86"/>
<instance part="U$15" gate="G$1" x="-78.74" y="-25.4"/>
<instance part="U$16" gate="G$1" x="-78.74" y="-27.94"/>
<instance part="U$17" gate="G$1" x="-78.74" y="-30.48"/>
<instance part="U$18" gate="G$1" x="-127" y="-30.48"/>
<instance part="U$19" gate="G$1" x="-127" y="-27.94"/>
<instance part="U$20" gate="G$1" x="-127" y="-25.4"/>
<instance part="U$21" gate="G$1" x="-127" y="-22.86"/>
<instance part="U$22" gate="G$1" x="-127" y="-20.32"/>
<instance part="U$23" gate="G$1" x="-127" y="-17.78"/>
<instance part="U$24" gate="G$1" x="-127" y="-15.24"/>
<instance part="U$25" gate="G$1" x="-127" y="-12.7"/>
<instance part="U$26" gate="G$1" x="-127" y="-10.16"/>
<instance part="U$27" gate="G$1" x="-127" y="-7.62"/>
<instance part="U$28" gate="G$1" x="-127" y="-5.08"/>
<instance part="U$29" gate="G$1" x="-127" y="-2.54"/>
<instance part="U$30" gate="G$1" x="-127" y="0"/>
<instance part="U$31" gate="G$1" x="-127" y="2.54"/>
<instance part="U$32" gate="G$1" x="-127" y="5.08"/>
<instance part="U$33" gate="G$1" x="-127" y="7.62"/>
<instance part="GND9" gate="1" x="60.96" y="81.28" rot="R270"/>
<instance part="SUPPLY1" gate="0V" x="60.96" y="76.2" rot="R270"/>
<instance part="U$34" gate="G$1" x="78.74" y="93.98" rot="R90"/>
<instance part="HEADSTAGE" gate="G$1" x="-5.08" y="-35.56" rot="MR270"/>
<instance part="U$35" gate="G$1" x="78.74" y="106.68" rot="R90"/>
<instance part="U$36" gate="G$1" x="78.74" y="119.38" rot="R90"/>
<instance part="SUPPLY2" gate="0V" x="30.48" y="101.6" rot="R90"/>
</instances>
<busses>
</busses>
<nets>
<net name="GND" class="0">
<segment>
<pinref part="GND1" gate="1" pin="GND"/>
<wire x1="-121.92" y1="-35.56" x2="-124.46" y2="-35.56" width="0.1524" layer="91"/>
<pinref part="STIM" gate="G$1" pin="35"/>
</segment>
<segment>
<pinref part="GND2" gate="1" pin="GND"/>
<wire x1="-86.36" y1="-33.02" x2="-81.28" y2="-33.02" width="0.1524" layer="91"/>
<pinref part="STIM" gate="G$1" pin="34"/>
</segment>
<segment>
<pinref part="GND3" gate="1" pin="GND"/>
<wire x1="17.78" y1="-17.78" x2="33.02" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="HEADSTAGE" gate="G$1" pin="35"/>
</segment>
<segment>
<wire x1="15.24" y1="-53.34" x2="30.48" y2="-58.42" width="0.1524" layer="91"/>
<pinref part="GND4" gate="1" pin="GND"/>
<wire x1="30.48" y1="-58.42" x2="33.02" y2="-58.42" width="0.1524" layer="91"/>
<pinref part="HEADSTAGE" gate="G$1" pin="34"/>
</segment>
<segment>
<pinref part="BUFFER" gate="G$1" pin="34"/>
<pinref part="GND5" gate="1" pin="GND"/>
<wire x1="2.54" y1="147.32" x2="5.08" y2="154.94" width="0.1524" layer="91"/>
</segment>
<segment>
<pinref part="BUFFER" gate="G$1" pin="35"/>
<pinref part="GND6" gate="1" pin="GND"/>
<wire x1="5.08" y1="111.76" x2="25.4" y2="111.76" width="0.1524" layer="91"/>
</segment>
<segment>
<pinref part="GND9" gate="1" pin="GND"/>
<wire x1="63.5" y1="81.28" x2="81.28" y2="81.28" width="0.1524" layer="91"/>
<pinref part="U$34" gate="G$1" pin="P$4"/>
<pinref part="U$35" gate="G$1" pin="P$4"/>
<wire x1="81.28" y1="88.9" x2="81.28" y2="101.6" width="0.1524" layer="91"/>
<pinref part="U$36" gate="G$1" pin="P$4"/>
<wire x1="81.28" y1="101.6" x2="81.28" y2="114.3" width="0.1524" layer="91"/>
<junction x="81.28" y="101.6"/>
<wire x1="81.28" y1="81.28" x2="81.28" y2="88.9" width="0.1524" layer="91"/>
<junction x="81.28" y="88.9"/>
</segment>
</net>
<net name="N$1" class="0">
<segment>
<wire x1="-38.1" y1="111.76" x2="-38.1" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="1"/>
<pinref part="U$3" gate="G$1" pin="2"/>
</segment>
</net>
<net name="N$3" class="0">
<segment>
<wire x1="-35.56" y1="111.76" x2="-33.02" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="3"/>
<pinref part="U$3" gate="G$1" pin="6"/>
</segment>
</net>
<net name="N$5" class="0">
<segment>
<wire x1="-33.02" y1="111.76" x2="-27.94" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="5"/>
<pinref part="U$3" gate="G$1" pin="10"/>
</segment>
</net>
<net name="N$6" class="0">
<segment>
<wire x1="-30.48" y1="111.76" x2="-22.86" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="7"/>
<pinref part="U$3" gate="G$1" pin="14"/>
</segment>
</net>
<net name="N$7" class="0">
<segment>
<wire x1="-27.94" y1="111.76" x2="-17.78" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="9"/>
<pinref part="U$3" gate="G$1" pin="18"/>
</segment>
</net>
<net name="N$8" class="0">
<segment>
<wire x1="-25.4" y1="111.76" x2="-12.7" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="11"/>
<pinref part="U$3" gate="G$1" pin="22"/>
</segment>
</net>
<net name="N$9" class="0">
<segment>
<wire x1="-7.62" y1="58.42" x2="-22.86" y2="111.76" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="13"/>
<pinref part="U$3" gate="G$1" pin="26"/>
</segment>
</net>
<net name="N$10" class="0">
<segment>
<wire x1="-20.32" y1="111.76" x2="-2.54" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="15"/>
<pinref part="U$3" gate="G$1" pin="30"/>
</segment>
</net>
<net name="N$11" class="0">
<segment>
<wire x1="2.54" y1="58.42" x2="-17.78" y2="111.76" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="17"/>
<pinref part="U$3" gate="G$1" pin="34"/>
</segment>
</net>
<net name="N$12" class="0">
<segment>
<wire x1="-15.24" y1="111.76" x2="7.62" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="19"/>
<pinref part="U$3" gate="G$1" pin="38"/>
</segment>
</net>
<net name="N$14" class="0">
<segment>
<wire x1="-10.16" y1="111.76" x2="17.78" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="23"/>
<pinref part="U$3" gate="G$1" pin="46"/>
</segment>
</net>
<net name="N$15" class="0">
<segment>
<wire x1="22.86" y1="58.42" x2="-7.62" y2="111.76" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="25"/>
<pinref part="U$3" gate="G$1" pin="50"/>
</segment>
</net>
<net name="N$18" class="0">
<segment>
<wire x1="0" y1="111.76" x2="38.1" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="31"/>
<pinref part="U$3" gate="G$1" pin="62"/>
</segment>
</net>
<net name="N$26" class="0">
<segment>
<wire x1="-12.7" y1="111.76" x2="12.7" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="21"/>
<pinref part="U$3" gate="G$1" pin="42"/>
</segment>
</net>
<net name="N$13" class="0">
<segment>
<wire x1="-5.08" y1="111.76" x2="27.94" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="27"/>
<pinref part="U$3" gate="G$1" pin="54"/>
</segment>
</net>
<net name="N$2" class="0">
<segment>
<wire x1="-2.54" y1="111.76" x2="33.02" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="29"/>
<pinref part="U$3" gate="G$1" pin="58"/>
</segment>
</net>
<net name="N$4" class="0">
<segment>
<wire x1="-35.56" y1="58.42" x2="-38.1" y2="147.32" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="2"/>
<pinref part="U$3" gate="G$1" pin="4"/>
</segment>
</net>
<net name="N$16" class="0">
<segment>
<wire x1="-35.56" y1="147.32" x2="-30.48" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="4"/>
<pinref part="U$3" gate="G$1" pin="8"/>
</segment>
</net>
<net name="N$17" class="0">
<segment>
<wire x1="-25.4" y1="58.42" x2="-33.02" y2="147.32" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="6"/>
<pinref part="U$3" gate="G$1" pin="12"/>
</segment>
</net>
<net name="N$19" class="0">
<segment>
<wire x1="-30.48" y1="147.32" x2="-20.32" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="8"/>
<pinref part="U$3" gate="G$1" pin="16"/>
</segment>
</net>
<net name="N$20" class="0">
<segment>
<wire x1="-27.94" y1="147.32" x2="-15.24" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="10"/>
<pinref part="U$3" gate="G$1" pin="20"/>
</segment>
</net>
<net name="N$21" class="0">
<segment>
<wire x1="-10.16" y1="58.42" x2="-25.4" y2="147.32" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="12"/>
<pinref part="U$3" gate="G$1" pin="24"/>
</segment>
</net>
<net name="N$22" class="0">
<segment>
<wire x1="-22.86" y1="147.32" x2="-5.08" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="14"/>
<pinref part="U$3" gate="G$1" pin="28"/>
</segment>
</net>
<net name="N$23" class="0">
<segment>
<wire x1="0" y1="58.42" x2="-20.32" y2="147.32" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="16"/>
<pinref part="U$3" gate="G$1" pin="32"/>
</segment>
</net>
<net name="N$24" class="0">
<segment>
<wire x1="-17.78" y1="147.32" x2="5.08" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="18"/>
<pinref part="U$3" gate="G$1" pin="36"/>
</segment>
</net>
<net name="N$25" class="0">
<segment>
<wire x1="10.16" y1="58.42" x2="-15.24" y2="147.32" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="20"/>
<pinref part="U$3" gate="G$1" pin="40"/>
</segment>
</net>
<net name="N$27" class="0">
<segment>
<wire x1="-12.7" y1="147.32" x2="15.24" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="22"/>
<pinref part="U$3" gate="G$1" pin="44"/>
</segment>
</net>
<net name="N$28" class="0">
<segment>
<wire x1="20.32" y1="58.42" x2="-10.16" y2="147.32" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="24"/>
<pinref part="U$3" gate="G$1" pin="48"/>
</segment>
</net>
<net name="N$29" class="0">
<segment>
<wire x1="-7.62" y1="147.32" x2="25.4" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="26"/>
<pinref part="U$3" gate="G$1" pin="52"/>
</segment>
</net>
<net name="N$30" class="0">
<segment>
<wire x1="30.48" y1="58.42" x2="-5.08" y2="147.32" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="28"/>
<pinref part="U$3" gate="G$1" pin="56"/>
</segment>
</net>
<net name="N$31" class="0">
<segment>
<wire x1="-2.54" y1="147.32" x2="35.56" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="30"/>
<pinref part="U$3" gate="G$1" pin="60"/>
</segment>
</net>
<net name="N$32" class="0">
<segment>
<wire x1="0" y1="147.32" x2="40.64" y2="58.42" width="0.1524" layer="91"/>
<pinref part="BUFFER" gate="G$1" pin="32"/>
<pinref part="U$3" gate="G$1" pin="64"/>
</segment>
</net>
<net name="N$33" class="0">
<segment>
<wire x1="-25.4" y1="-17.78" x2="-38.1" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="1"/>
<pinref part="HEADSTAGE" gate="G$1" pin="1"/>
</segment>
</net>
<net name="N$34" class="0">
<segment>
<wire x1="-33.02" y1="22.86" x2="-22.86" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="5"/>
<pinref part="HEADSTAGE" gate="G$1" pin="3"/>
</segment>
</net>
<net name="N$35" class="0">
<segment>
<wire x1="-20.32" y1="-17.78" x2="-27.94" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="9"/>
<pinref part="HEADSTAGE" gate="G$1" pin="5"/>
</segment>
</net>
<net name="N$36" class="0">
<segment>
<wire x1="-22.86" y1="22.86" x2="-17.78" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="13"/>
<pinref part="HEADSTAGE" gate="G$1" pin="7"/>
</segment>
</net>
<net name="N$37" class="0">
<segment>
<wire x1="-15.24" y1="-17.78" x2="-17.78" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="17"/>
<pinref part="HEADSTAGE" gate="G$1" pin="9"/>
</segment>
</net>
<net name="N$38" class="0">
<segment>
<wire x1="-12.7" y1="22.86" x2="-12.7" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="21"/>
<pinref part="HEADSTAGE" gate="G$1" pin="11"/>
</segment>
</net>
<net name="N$39" class="0">
<segment>
<wire x1="-10.16" y1="-17.78" x2="-7.62" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="25"/>
<pinref part="HEADSTAGE" gate="G$1" pin="13"/>
</segment>
</net>
<net name="N$40" class="0">
<segment>
<wire x1="-2.54" y1="22.86" x2="-7.62" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="29"/>
<pinref part="HEADSTAGE" gate="G$1" pin="15"/>
</segment>
</net>
<net name="N$41" class="0">
<segment>
<wire x1="-5.08" y1="-17.78" x2="2.54" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="33"/>
<pinref part="HEADSTAGE" gate="G$1" pin="17"/>
</segment>
</net>
<net name="N$42" class="0">
<segment>
<wire x1="7.62" y1="22.86" x2="-2.54" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="37"/>
<pinref part="HEADSTAGE" gate="G$1" pin="19"/>
</segment>
</net>
<net name="N$43" class="0">
<segment>
<wire x1="0" y1="-17.78" x2="12.7" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="41"/>
<pinref part="HEADSTAGE" gate="G$1" pin="21"/>
</segment>
</net>
<net name="N$44" class="0">
<segment>
<wire x1="17.78" y1="22.86" x2="2.54" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="45"/>
<pinref part="HEADSTAGE" gate="G$1" pin="23"/>
</segment>
</net>
<net name="N$45" class="0">
<segment>
<wire x1="5.08" y1="-17.78" x2="22.86" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="49"/>
<pinref part="HEADSTAGE" gate="G$1" pin="25"/>
</segment>
</net>
<net name="N$46" class="0">
<segment>
<wire x1="27.94" y1="22.86" x2="7.62" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="53"/>
<pinref part="HEADSTAGE" gate="G$1" pin="27"/>
</segment>
</net>
<net name="N$47" class="0">
<segment>
<wire x1="10.16" y1="-17.78" x2="33.02" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="57"/>
<pinref part="HEADSTAGE" gate="G$1" pin="29"/>
</segment>
</net>
<net name="N$48" class="0">
<segment>
<wire x1="38.1" y1="22.86" x2="12.7" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="61"/>
<pinref part="HEADSTAGE" gate="G$1" pin="31"/>
</segment>
</net>
<net name="N$49" class="0">
<segment>
<wire x1="-35.56" y1="22.86" x2="-25.4" y2="-53.34" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="3"/>
<pinref part="HEADSTAGE" gate="G$1" pin="2"/>
</segment>
</net>
<net name="N$50" class="0">
<segment>
<wire x1="-22.86" y1="-53.34" x2="-30.48" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="7"/>
<pinref part="HEADSTAGE" gate="G$1" pin="4"/>
</segment>
</net>
<net name="N$51" class="0">
<segment>
<wire x1="-25.4" y1="22.86" x2="-20.32" y2="-53.34" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="11"/>
<pinref part="HEADSTAGE" gate="G$1" pin="6"/>
</segment>
</net>
<net name="N$52" class="0">
<segment>
<wire x1="-17.78" y1="-53.34" x2="-20.32" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="15"/>
<pinref part="HEADSTAGE" gate="G$1" pin="8"/>
</segment>
</net>
<net name="N$53" class="0">
<segment>
<wire x1="-15.24" y1="22.86" x2="-15.24" y2="-53.34" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="19"/>
<pinref part="HEADSTAGE" gate="G$1" pin="10"/>
</segment>
</net>
<net name="N$54" class="0">
<segment>
<wire x1="-12.7" y1="-53.34" x2="-10.16" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="23"/>
<pinref part="HEADSTAGE" gate="G$1" pin="12"/>
</segment>
</net>
<net name="N$55" class="0">
<segment>
<wire x1="-5.08" y1="22.86" x2="-10.16" y2="-53.34" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="27"/>
<pinref part="HEADSTAGE" gate="G$1" pin="14"/>
</segment>
</net>
<net name="N$56" class="0">
<segment>
<wire x1="-7.62" y1="-53.34" x2="0" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="31"/>
<pinref part="HEADSTAGE" gate="G$1" pin="16"/>
</segment>
</net>
<net name="N$57" class="0">
<segment>
<wire x1="5.08" y1="22.86" x2="-5.08" y2="-53.34" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="35"/>
<pinref part="HEADSTAGE" gate="G$1" pin="18"/>
</segment>
</net>
<net name="N$58" class="0">
<segment>
<wire x1="-2.54" y1="-53.34" x2="10.16" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="39"/>
<pinref part="HEADSTAGE" gate="G$1" pin="20"/>
</segment>
</net>
<net name="N$59" class="0">
<segment>
<wire x1="15.24" y1="22.86" x2="0" y2="-53.34" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="43"/>
<pinref part="HEADSTAGE" gate="G$1" pin="22"/>
</segment>
</net>
<net name="N$60" class="0">
<segment>
<wire x1="2.54" y1="-53.34" x2="20.32" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="47"/>
<pinref part="HEADSTAGE" gate="G$1" pin="24"/>
</segment>
</net>
<net name="N$61" class="0">
<segment>
<wire x1="25.4" y1="22.86" x2="5.08" y2="-53.34" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="51"/>
<pinref part="HEADSTAGE" gate="G$1" pin="26"/>
</segment>
</net>
<net name="N$62" class="0">
<segment>
<wire x1="7.62" y1="-53.34" x2="30.48" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="55"/>
<pinref part="HEADSTAGE" gate="G$1" pin="28"/>
</segment>
</net>
<net name="N$63" class="0">
<segment>
<wire x1="35.56" y1="22.86" x2="10.16" y2="-53.34" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="59"/>
<pinref part="HEADSTAGE" gate="G$1" pin="30"/>
</segment>
</net>
<net name="N$64" class="0">
<segment>
<wire x1="12.7" y1="-53.34" x2="40.64" y2="22.86" width="0.1524" layer="91"/>
<pinref part="U$3" gate="G$1" pin="63"/>
<pinref part="HEADSTAGE" gate="G$1" pin="32"/>
</segment>
</net>
<net name="0V" class="0">
<segment>
<pinref part="GND7" gate="1" pin="0V"/>
<wire x1="33.02" y1="-10.16" x2="15.24" y2="-17.78" width="0.1524" layer="91"/>
<pinref part="HEADSTAGE" gate="G$1" pin="33"/>
</segment>
<segment>
<pinref part="SUPPLY1" gate="0V" pin="0V"/>
<wire x1="63.5" y1="76.2" x2="76.2" y2="88.9" width="0.1524" layer="91"/>
<pinref part="U$34" gate="G$1" pin="P$2"/>
<pinref part="U$35" gate="G$1" pin="P$2"/>
<wire x1="76.2" y1="88.9" x2="76.2" y2="101.6" width="0.1524" layer="91"/>
<junction x="76.2" y="88.9"/>
<pinref part="U$36" gate="G$1" pin="P$2"/>
<wire x1="76.2" y1="101.6" x2="76.2" y2="114.3" width="0.1524" layer="91"/>
<junction x="76.2" y="101.6"/>
</segment>
<segment>
<pinref part="BUFFER" gate="G$1" pin="33"/>
<pinref part="SUPPLY2" gate="0V" pin="0V"/>
<wire x1="2.54" y1="111.76" x2="27.94" y2="101.6" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$65" class="0">
<segment>
<pinref part="U$33" gate="G$1" pin="P$1"/>
<pinref part="STIM" gate="G$1" pin="1"/>
<wire x1="-132.08" y1="7.62" x2="-121.92" y2="7.62" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$66" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="3"/>
<pinref part="U$32" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="5.08" x2="-132.08" y2="5.08" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$67" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="5"/>
<pinref part="U$31" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="2.54" x2="-132.08" y2="2.54" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$68" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="7"/>
<pinref part="U$30" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="0" x2="-132.08" y2="0" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$69" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="9"/>
<pinref part="U$29" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-2.54" x2="-132.08" y2="-2.54" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$70" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="11"/>
<pinref part="U$28" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-5.08" x2="-132.08" y2="-5.08" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$71" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="13"/>
<pinref part="U$27" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-7.62" x2="-132.08" y2="-7.62" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$72" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="15"/>
<pinref part="U$26" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-10.16" x2="-132.08" y2="-10.16" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$73" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="17"/>
<pinref part="U$25" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-12.7" x2="-132.08" y2="-12.7" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$74" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="19"/>
<pinref part="U$24" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-15.24" x2="-132.08" y2="-15.24" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$75" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="21"/>
<pinref part="U$23" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-17.78" x2="-132.08" y2="-17.78" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$76" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="23"/>
<pinref part="U$22" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-20.32" x2="-132.08" y2="-20.32" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$77" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="25"/>
<pinref part="U$21" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-22.86" x2="-132.08" y2="-22.86" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$78" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="27"/>
<pinref part="U$20" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-25.4" x2="-132.08" y2="-25.4" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$79" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="29"/>
<pinref part="U$19" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-27.94" x2="-132.08" y2="-27.94" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$80" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="31"/>
<pinref part="U$18" gate="G$1" pin="P$1"/>
<wire x1="-121.92" y1="-30.48" x2="-132.08" y2="-30.48" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$81" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="32"/>
<pinref part="U$17" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-30.48" x2="-83.82" y2="-30.48" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$82" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="30"/>
<pinref part="U$16" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-27.94" x2="-83.82" y2="-27.94" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$83" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="28"/>
<pinref part="U$15" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-25.4" x2="-83.82" y2="-25.4" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$84" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="26"/>
<pinref part="U$14" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-22.86" x2="-83.82" y2="-22.86" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$85" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="24"/>
<pinref part="U$13" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-20.32" x2="-83.82" y2="-20.32" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$86" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="22"/>
<pinref part="U$12" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-17.78" x2="-83.82" y2="-17.78" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$87" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="20"/>
<pinref part="U$11" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-15.24" x2="-83.82" y2="-15.24" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$88" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="18"/>
<pinref part="U$10" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-12.7" x2="-83.82" y2="-12.7" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$89" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="16"/>
<pinref part="U$9" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-10.16" x2="-83.82" y2="-10.16" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$90" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="14"/>
<pinref part="U$8" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-7.62" x2="-83.82" y2="-7.62" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$91" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="12"/>
<pinref part="U$7" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-5.08" x2="-83.82" y2="-5.08" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$92" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="10"/>
<pinref part="U$6" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="-2.54" x2="-83.82" y2="-2.54" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$93" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="8"/>
<pinref part="U$5" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="0" x2="-83.82" y2="0" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$94" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="6"/>
<pinref part="U$4" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="2.54" x2="-83.82" y2="2.54" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$95" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="4"/>
<pinref part="U$2" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="5.08" x2="-83.82" y2="5.08" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$96" class="0">
<segment>
<pinref part="STIM" gate="G$1" pin="2"/>
<pinref part="U$1" gate="G$1" pin="P$1"/>
<wire x1="-86.36" y1="7.62" x2="-83.82" y2="7.62" width="0.1524" layer="91"/>
</segment>
</net>
<net name="N$97" class="0">
<segment>
<wire x1="96.52" y1="76.2" x2="96.52" y2="58.42" width="0.1524" layer="91"/>
<pinref part="DUKE-STIM" gate="-1" pin="S"/>
<wire x1="96.52" y1="58.42" x2="101.6" y2="58.42" width="0.1524" layer="91"/>
<pinref part="DUKE-REC" gate="-1" pin="S"/>
<pinref part="U$34" gate="G$1" pin="P$3"/>
<wire x1="78.74" y1="88.9" x2="78.74" y2="76.2" width="0.1524" layer="91"/>
<wire x1="78.74" y1="76.2" x2="96.52" y2="76.2" width="0.1524" layer="91"/>
<wire x1="96.52" y1="76.2" x2="99.06" y2="76.2" width="0.1524" layer="91"/>
<junction x="96.52" y="76.2"/>
<pinref part="U$35" gate="G$1" pin="P$3"/>
<wire x1="78.74" y1="88.9" x2="78.74" y2="101.6" width="0.1524" layer="91"/>
<junction x="78.74" y="88.9"/>
<pinref part="U$36" gate="G$1" pin="P$3"/>
<wire x1="78.74" y1="101.6" x2="78.74" y2="114.3" width="0.1524" layer="91"/>
<junction x="78.74" y="101.6"/>
</segment>
</net>
<net name="N$99" class="0">
<segment>
<wire x1="83.82" y1="88.9" x2="83.82" y2="55.88" width="0.1524" layer="91"/>
<pinref part="DUKE-STIM" gate="-2" pin="S"/>
<wire x1="83.82" y1="55.88" x2="101.6" y2="55.88" width="0.1524" layer="91"/>
<pinref part="U$34" gate="G$1" pin="P$5"/>
<pinref part="U$35" gate="G$1" pin="P$5"/>
<wire x1="83.82" y1="88.9" x2="83.82" y2="101.6" width="0.1524" layer="91"/>
<junction x="83.82" y="88.9"/>
<pinref part="U$36" gate="G$1" pin="P$5"/>
<wire x1="83.82" y1="101.6" x2="83.82" y2="114.3" width="0.1524" layer="91"/>
<junction x="83.82" y="101.6"/>
</segment>
</net>
<net name="N$98" class="0">
<segment>
<pinref part="U$34" gate="G$1" pin="P$1"/>
<pinref part="U$35" gate="G$1" pin="P$1"/>
<wire x1="73.66" y1="88.9" x2="73.66" y2="101.6" width="0.1524" layer="91"/>
<pinref part="U$36" gate="G$1" pin="P$1"/>
<wire x1="73.66" y1="101.6" x2="73.66" y2="114.3" width="0.1524" layer="91"/>
<junction x="73.66" y="101.6"/>
<wire x1="73.66" y1="88.9" x2="73.66" y2="73.66" width="0.1524" layer="91"/>
<junction x="73.66" y="88.9"/>
<pinref part="DUKE-REC" gate="-2" pin="S"/>
<wire x1="73.66" y1="73.66" x2="99.06" y2="73.66" width="0.1524" layer="91"/>
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
