% Calculate Gravitation load for parallel robot
% P3PRRRR8V1G4A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:24:57
% EndTime: 2020-08-06 17:24:59
% DurationCPUTime: 2.01s
% Computational Cost: add. (1623->206), mult. (4128->426), div. (54->7), fcn. (4263->34), ass. (0->168)
t1164 = sin(qJ(2,3));
t1170 = cos(qJ(2,3));
t1154 = legFrame(3,3);
t1128 = sin(t1154);
t1157 = legFrame(3,1);
t1131 = sin(t1157);
t1134 = cos(t1154);
t1137 = cos(t1157);
t1160 = legFrame(3,2);
t1140 = sin(t1160);
t1218 = t1137 * t1140;
t1221 = t1131 * t1140;
t1143 = cos(t1160);
t1233 = g(1) * t1143;
t1061 = -t1128 * t1233 + (-t1128 * t1221 + t1134 * t1137) * g(2) + (t1128 * t1218 + t1131 * t1134) * g(3);
t1062 = t1134 * t1233 + (t1128 * t1137 + t1134 * t1221) * g(2) + (t1128 * t1131 - t1134 * t1218) * g(3);
t1150 = sin(pkin(6));
t1152 = cos(pkin(6));
t1193 = t1061 * t1150 + t1062 * t1152;
t1088 = g(1) * t1140 + (-g(2) * t1131 + g(3) * t1137) * t1143;
t1151 = sin(pkin(3));
t1153 = cos(pkin(3));
t1194 = t1061 * t1152 - t1062 * t1150;
t1244 = t1088 * t1151 + t1194 * t1153;
t1179 = t1244 * t1164 + t1193 * t1170;
t1166 = sin(qJ(2,2));
t1172 = cos(qJ(2,2));
t1155 = legFrame(2,3);
t1129 = sin(t1155);
t1158 = legFrame(2,1);
t1132 = sin(t1158);
t1135 = cos(t1155);
t1138 = cos(t1158);
t1161 = legFrame(2,2);
t1141 = sin(t1161);
t1217 = t1138 * t1141;
t1220 = t1132 * t1141;
t1144 = cos(t1161);
t1232 = g(1) * t1144;
t1063 = -t1129 * t1232 + (-t1129 * t1220 + t1135 * t1138) * g(2) + (t1129 * t1217 + t1132 * t1135) * g(3);
t1064 = t1135 * t1232 + (t1129 * t1138 + t1135 * t1220) * g(2) + (t1129 * t1132 - t1135 * t1217) * g(3);
t1191 = t1063 * t1150 + t1064 * t1152;
t1089 = g(1) * t1141 + (-g(2) * t1132 + g(3) * t1138) * t1144;
t1192 = t1063 * t1152 - t1064 * t1150;
t1245 = t1089 * t1151 + t1192 * t1153;
t1178 = t1245 * t1166 + t1191 * t1172;
t1168 = sin(qJ(2,1));
t1174 = cos(qJ(2,1));
t1156 = legFrame(1,3);
t1130 = sin(t1156);
t1159 = legFrame(1,1);
t1133 = sin(t1159);
t1136 = cos(t1156);
t1139 = cos(t1159);
t1162 = legFrame(1,2);
t1142 = sin(t1162);
t1216 = t1139 * t1142;
t1219 = t1133 * t1142;
t1145 = cos(t1162);
t1231 = g(1) * t1145;
t1065 = -t1130 * t1231 + (-t1130 * t1219 + t1136 * t1139) * g(2) + (t1130 * t1216 + t1133 * t1136) * g(3);
t1066 = t1136 * t1231 + (t1130 * t1139 + t1136 * t1219) * g(2) + (t1130 * t1133 - t1136 * t1216) * g(3);
t1189 = t1065 * t1150 + t1066 * t1152;
t1090 = g(1) * t1142 + (-g(2) * t1133 + g(3) * t1139) * t1145;
t1190 = t1065 * t1152 - t1066 * t1150;
t1246 = t1090 * t1151 + t1190 * t1153;
t1177 = t1246 * t1168 + t1189 * t1174;
t1169 = cos(qJ(3,3));
t1236 = pkin(2) * t1169;
t1115 = -pkin(5) * t1170 + t1164 * t1236;
t1163 = sin(qJ(3,3));
t1239 = pkin(2) * t1163;
t1097 = t1115 * t1151 + t1153 * t1239;
t1243 = 0.1e1 / t1097;
t1252 = 0.1e1 / t1169 * t1243;
t1171 = cos(qJ(3,2));
t1235 = pkin(2) * t1171;
t1116 = -pkin(5) * t1172 + t1166 * t1235;
t1165 = sin(qJ(3,2));
t1238 = pkin(2) * t1165;
t1098 = t1116 * t1151 + t1153 * t1238;
t1242 = 0.1e1 / t1098;
t1251 = 0.1e1 / t1171 * t1242;
t1173 = cos(qJ(3,1));
t1234 = pkin(2) * t1173;
t1117 = -pkin(5) * t1174 + t1168 * t1234;
t1167 = sin(qJ(3,1));
t1237 = pkin(2) * t1167;
t1099 = t1117 * t1151 + t1153 * t1237;
t1241 = 0.1e1 / t1099;
t1250 = 0.1e1 / t1173 * t1241;
t1249 = -t1090 * t1153 + t1190 * t1151;
t1248 = -t1089 * t1153 + t1192 * t1151;
t1247 = -t1088 * t1153 + t1194 * t1151;
t1240 = m(3) / pkin(2);
t1230 = t1088 * t1243;
t1227 = t1089 * t1242;
t1224 = t1090 * t1241;
t1209 = t1153 * t1164;
t1106 = t1150 * t1170 + t1152 * t1209;
t1109 = -t1150 * t1209 + t1152 * t1170;
t1085 = t1106 * t1134 + t1109 * t1128;
t1215 = t1140 * t1085;
t1208 = t1153 * t1166;
t1107 = t1150 * t1172 + t1152 * t1208;
t1110 = -t1150 * t1208 + t1152 * t1172;
t1086 = t1107 * t1135 + t1110 * t1129;
t1214 = t1141 * t1086;
t1207 = t1153 * t1168;
t1108 = t1150 * t1174 + t1152 * t1207;
t1111 = -t1150 * t1207 + t1152 * t1174;
t1087 = t1108 * t1136 + t1111 * t1130;
t1213 = t1142 * t1087;
t1212 = t1151 * t1169;
t1211 = t1151 * t1171;
t1210 = t1151 * t1173;
t1206 = t1153 * t1170;
t1205 = t1153 * t1172;
t1204 = t1153 * t1174;
t1203 = ((t1247 * rSges(3,1) + t1179 * rSges(3,2)) * t1169 + t1163 * (t1179 * rSges(3,1) - t1247 * rSges(3,2))) * t1252;
t1202 = ((t1248 * rSges(3,1) + t1178 * rSges(3,2)) * t1171 + t1165 * (t1178 * rSges(3,1) - t1248 * rSges(3,2))) * t1251;
t1201 = ((t1249 * rSges(3,1) + t1177 * rSges(3,2)) * t1173 + t1167 * (t1177 * rSges(3,1) - t1249 * rSges(3,2))) * t1250;
t1127 = m(2) * rSges(2,2) - m(3) * rSges(3,3);
t1175 = m(2) * rSges(2,1);
t1200 = (t1179 * t1127 + (t1164 * t1193 - t1170 * t1244) * (t1175 + (rSges(3,1) * t1169 - rSges(3,2) * t1163) * m(3))) * t1252;
t1199 = (t1178 * t1127 + (t1166 * t1191 - t1172 * t1245) * (t1175 + (rSges(3,1) * t1171 - rSges(3,2) * t1165) * m(3))) * t1251;
t1198 = (t1177 * t1127 + (t1168 * t1189 - t1174 * t1246) * (t1175 + (rSges(3,1) * t1173 - rSges(3,2) * t1167) * m(3))) * t1250;
t1188 = t1106 * t1128 - t1109 * t1134;
t1187 = t1107 * t1129 - t1110 * t1135;
t1186 = t1108 * t1130 - t1111 * t1136;
t1185 = -t1115 * t1153 + t1151 * t1239;
t1184 = -t1116 * t1153 + t1151 * t1238;
t1183 = -t1117 * t1153 + t1151 * t1237;
t1146 = m(1) + m(2) + m(3);
t1120 = pkin(5) * t1168 + t1174 * t1234;
t1119 = pkin(5) * t1166 + t1172 * t1235;
t1118 = pkin(5) * t1164 + t1170 * t1236;
t1105 = -t1130 * t1150 + t1136 * t1152;
t1104 = -t1129 * t1150 + t1135 * t1152;
t1103 = -t1128 * t1150 + t1134 * t1152;
t1102 = t1130 * t1152 + t1136 * t1150;
t1101 = t1129 * t1152 + t1135 * t1150;
t1100 = t1128 * t1152 + t1134 * t1150;
t1084 = -t1120 * t1150 + t1183 * t1152;
t1083 = -t1119 * t1150 + t1184 * t1152;
t1082 = -t1118 * t1150 + t1185 * t1152;
t1081 = t1120 * t1152 + t1183 * t1150;
t1080 = t1119 * t1152 + t1184 * t1150;
t1079 = t1118 * t1152 + t1185 * t1150;
t1078 = t1102 * t1219 - t1105 * t1139;
t1077 = t1101 * t1220 - t1104 * t1138;
t1076 = t1100 * t1221 - t1103 * t1137;
t1075 = t1102 * t1139 + t1105 * t1219;
t1074 = t1101 * t1138 + t1104 * t1220;
t1073 = t1100 * t1137 + t1103 * t1221;
t1072 = -t1102 * t1133 + t1105 * t1216;
t1071 = t1102 * t1216 + t1105 * t1133;
t1070 = -t1101 * t1132 + t1104 * t1217;
t1069 = t1101 * t1217 + t1104 * t1132;
t1068 = -t1100 * t1131 + t1103 * t1218;
t1067 = t1100 * t1218 + t1103 * t1131;
t1060 = -t1081 * t1130 + t1084 * t1136;
t1059 = -t1080 * t1129 + t1083 * t1135;
t1058 = -t1079 * t1128 + t1082 * t1134;
t1054 = -t1099 * t1145 + (t1081 * t1136 + t1084 * t1130) * t1142;
t1053 = -t1098 * t1144 + (t1080 * t1135 + t1083 * t1129) * t1141;
t1052 = -t1097 * t1143 + (t1079 * t1134 + t1082 * t1128) * t1140;
t1 = [-(t1087 * t1167 + t1105 * t1210) * t1145 * t1198 - (t1086 * t1165 + t1104 * t1211) * t1144 * t1199 - (t1085 * t1163 + t1103 * t1212) * t1143 * t1200 - m(4) * g(1) + (-((t1183 * t1102 + t1105 * t1120) * t1145 + t1142 * t1099) * t1224 - ((t1184 * t1101 + t1104 * t1119) * t1144 + t1141 * t1098) * t1227 - ((t1185 * t1100 + t1103 * t1118) * t1143 + t1140 * t1097) * t1230) * t1146 + (-t1145 * ((-t1102 * t1168 + t1105 * t1204) * t1234 + pkin(5) * (t1102 * t1174 + t1105 * t1207)) * t1201 - t1144 * ((-t1101 * t1166 + t1104 * t1205) * t1235 + pkin(5) * (t1101 * t1172 + t1104 * t1208)) * t1202 - t1143 * ((-t1100 * t1164 + t1103 * t1206) * t1236 + pkin(5) * (t1100 * t1170 + t1103 * t1209)) * t1203) * t1240; ((-t1133 * t1213 - t1186 * t1139) * t1167 - t1075 * t1210) * t1198 + ((-t1132 * t1214 - t1187 * t1138) * t1165 - t1074 * t1211) * t1199 + ((-t1131 * t1215 - t1188 * t1137) * t1163 - t1073 * t1212) * t1200 - m(4) * g(2) + (-(t1054 * t1133 - t1060 * t1139) * t1224 - (t1053 * t1132 - t1059 * t1138) * t1227 - (t1052 * t1131 - t1058 * t1137) * t1230) * t1146 + ((-(t1075 * t1204 - t1078 * t1168) * t1234 - pkin(5) * (t1075 * t1207 + t1078 * t1174)) * t1201 + (-(t1074 * t1205 - t1077 * t1166) * t1235 - pkin(5) * (t1074 * t1208 + t1077 * t1172)) * t1202 + (-(t1073 * t1206 - t1076 * t1164) * t1236 - pkin(5) * (t1073 * t1209 + t1076 * t1170)) * t1203) * t1240; ((-t1186 * t1133 + t1139 * t1213) * t1167 + t1072 * t1210) * t1198 + ((-t1187 * t1132 + t1138 * t1214) * t1165 + t1070 * t1211) * t1199 + ((-t1188 * t1131 + t1137 * t1215) * t1163 + t1068 * t1212) * t1200 - m(4) * g(3) + (-(-t1054 * t1139 - t1060 * t1133) * t1224 - (-t1053 * t1138 - t1059 * t1132) * t1227 - (-t1052 * t1137 - t1058 * t1131) * t1230) * t1146 + (((-t1071 * t1168 + t1072 * t1204) * t1234 + pkin(5) * (t1071 * t1174 + t1072 * t1207)) * t1201 + ((-t1069 * t1166 + t1070 * t1205) * t1235 + pkin(5) * (t1069 * t1172 + t1070 * t1208)) * t1202 + ((-t1067 * t1164 + t1068 * t1206) * t1236 + pkin(5) * (t1067 * t1170 + t1068 * t1209)) * t1203) * t1240;];
taugX  = t1;
