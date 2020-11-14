% Calculate Gravitation load for parallel robot
% P3RRPRR8V2G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:11:50
% EndTime: 2020-08-06 21:11:51
% DurationCPUTime: 1.20s
% Computational Cost: add. (789->173), mult. (1449->282), div. (36->9), fcn. (792->23), ass. (0->128)
t1259 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1160 = sin(pkin(7));
t1161 = cos(pkin(7));
t1109 = -m(2) * rSges(2,1) + (-rSges(3,1) * t1161 + rSges(3,2) * t1160 - pkin(2)) * m(3);
t1111 = m(2) * rSges(2,2) + (rSges(3,1) * t1160 + rSges(3,2) * t1161) * m(3);
t1169 = sin(qJ(2,1));
t1175 = cos(qJ(2,1));
t1257 = -t1109 * t1175 - t1111 * t1169;
t1167 = sin(qJ(2,2));
t1173 = cos(qJ(2,2));
t1256 = -t1109 * t1173 - t1111 * t1167;
t1165 = sin(qJ(2,3));
t1171 = cos(qJ(2,3));
t1255 = -t1109 * t1171 - t1111 * t1165;
t1254 = 0.2e1 * t1171 ^ 2;
t1253 = 0.2e1 * t1173 ^ 2;
t1252 = 0.2e1 * t1175 ^ 2;
t1251 = pkin(5) + qJ(3,1);
t1250 = pkin(5) + qJ(3,2);
t1249 = pkin(5) + qJ(3,3);
t1166 = sin(qJ(1,3));
t1153 = pkin(6) + t1249;
t1172 = cos(qJ(1,3));
t1188 = pkin(1) * t1166 - t1172 * t1153;
t1156 = t1161 ^ 2;
t1209 = pkin(3) * (t1156 - 0.1e1);
t1213 = t1160 * t1165;
t1248 = pkin(3) * (t1166 * t1209 + t1188 * t1213);
t1168 = sin(qJ(1,2));
t1154 = pkin(6) + t1250;
t1174 = cos(qJ(1,2));
t1187 = pkin(1) * t1168 - t1174 * t1154;
t1212 = t1160 * t1167;
t1247 = pkin(3) * (t1168 * t1209 + t1187 * t1212);
t1170 = sin(qJ(1,1));
t1155 = pkin(6) + t1251;
t1176 = cos(qJ(1,1));
t1186 = pkin(1) * t1170 - t1176 * t1155;
t1211 = t1160 * t1169;
t1246 = pkin(3) * (t1170 * t1209 + t1186 * t1211);
t1245 = (rSges(3,3) + t1249) * m(3);
t1244 = (rSges(3,3) + t1250) * m(3);
t1243 = (rSges(3,3) + t1251) * m(3);
t1242 = t1160 * pkin(3);
t1241 = t1161 * pkin(3);
t1162 = legFrame(3,2);
t1137 = sin(t1162);
t1140 = cos(t1162);
t1115 = t1140 * g(1) - t1137 * g(2);
t1134 = 0.1e1 / t1153;
t1240 = (g(3) * t1166 - t1115 * t1172) * t1134;
t1163 = legFrame(2,2);
t1138 = sin(t1163);
t1141 = cos(t1163);
t1116 = t1141 * g(1) - t1138 * g(2);
t1135 = 0.1e1 / t1154;
t1239 = (g(3) * t1168 - t1116 * t1174) * t1135;
t1164 = legFrame(1,2);
t1139 = sin(t1164);
t1142 = cos(t1164);
t1117 = t1142 * g(1) - t1139 * g(2);
t1136 = 0.1e1 / t1155;
t1238 = (g(3) * t1170 - t1117 * t1176) * t1136;
t1112 = t1137 * g(1) + t1140 * g(2);
t1185 = g(3) * t1172 + t1115 * t1166;
t1191 = pkin(3) * cos(qJ(2,3) + pkin(7)) + t1171 * pkin(2);
t1231 = 0.1e1 / t1191 * ((-t1185 * t1109 + t1112 * t1111) * t1165 + (t1112 * t1109 + t1185 * t1111) * t1171);
t1113 = t1138 * g(1) + t1141 * g(2);
t1184 = g(3) * t1174 + t1116 * t1168;
t1190 = pkin(3) * cos(qJ(2,2) + pkin(7)) + t1173 * pkin(2);
t1230 = 0.1e1 / t1190 * ((-t1184 * t1109 + t1113 * t1111) * t1167 + (t1113 * t1109 + t1184 * t1111) * t1173);
t1114 = t1139 * g(1) + t1142 * g(2);
t1183 = g(3) * t1176 + t1117 * t1170;
t1189 = pkin(3) * cos(qJ(2,1) + pkin(7)) + t1175 * pkin(2);
t1229 = 0.1e1 / t1189 * ((-t1183 * t1109 + t1114 * t1111) * t1169 + (t1114 * t1109 + t1183 * t1111) * t1175);
t1127 = pkin(2) + t1241;
t1228 = t1127 * t1140;
t1227 = t1127 * t1141;
t1226 = t1127 * t1142;
t1125 = (-rSges(2,3) - pkin(5)) * m(2) + m(1) * rSges(1,2);
t1124 = t1125 * g(3);
t1192 = t1259 * g(3);
t1225 = t1134 * ((-g(3) * t1245 + t1124) * t1172 + t1166 * (t1255 * g(3) + t1192) + ((-t1259 - t1255) * t1172 + t1166 * (t1125 - t1245)) * t1115);
t1224 = t1135 * ((-g(3) * t1244 + t1124) * t1174 + t1168 * (t1256 * g(3) + t1192) + ((-t1259 - t1256) * t1174 + t1168 * (t1125 - t1244)) * t1116);
t1223 = t1136 * ((-g(3) * t1243 + t1124) * t1176 + t1170 * (t1257 * g(3) + t1192) + ((-t1259 - t1257) * t1176 + t1170 * (t1125 - t1243)) * t1117);
t1222 = t1137 * t1127;
t1221 = t1137 * t1166;
t1220 = t1138 * t1127;
t1219 = t1138 * t1168;
t1218 = t1139 * t1127;
t1217 = t1139 * t1170;
t1216 = t1140 * t1166;
t1215 = t1141 * t1168;
t1214 = t1142 * t1170;
t1210 = pkin(2) * t1241;
t1208 = t1140 * t1242;
t1207 = t1141 * t1242;
t1206 = t1142 * t1242;
t1205 = t1137 * t1242;
t1204 = t1138 * t1242;
t1203 = t1139 * t1242;
t1202 = pkin(3) * t1213;
t1201 = pkin(3) * t1212;
t1200 = pkin(3) * t1211;
t1106 = 0.1e1 / (t1127 * t1171 - t1202);
t1198 = t1106 * t1225;
t1107 = 0.1e1 / (t1127 * t1173 - t1201);
t1197 = t1107 * t1224;
t1108 = 0.1e1 / (t1127 * t1175 - t1200);
t1196 = t1108 * t1223;
t1195 = t1106 * t1240;
t1194 = t1107 * t1239;
t1193 = t1108 * t1238;
t1180 = pkin(3) ^ 2;
t1181 = pkin(2) ^ 2;
t1182 = 0.2e1 * t1156 * t1180 - t1180 + t1181 + 0.2e1 * t1210;
t1130 = pkin(1) * t1242;
t1123 = pkin(1) * t1169 - t1242;
t1122 = pkin(1) * t1167 - t1242;
t1121 = pkin(1) * t1165 - t1242;
t1110 = t1210 + t1181 / 0.2e1 + (t1156 - 0.1e1 / 0.2e1) * t1180;
t1102 = -0.2e1 * t1170 * t1200 + t1186;
t1101 = -0.2e1 * t1168 * t1201 + t1187;
t1100 = -0.2e1 * t1166 * t1202 + t1188;
t1099 = t1182 * t1169 + t1130;
t1098 = t1182 * t1167 + t1130;
t1097 = t1182 * t1165 + t1130;
t1 = [((t1127 * t1214 + t1203) * t1175 + t1169 * (-t1170 * t1206 + t1218)) * t1196 + t1139 * t1229 + ((t1127 * t1215 + t1204) * t1173 + t1167 * (-t1168 * t1207 + t1220)) * t1197 + t1138 * t1230 + ((t1127 * t1216 + t1205) * t1171 + t1165 * (-t1166 * t1208 + t1222)) * t1198 + t1137 * t1231 - m(4) * g(1) + (-((t1110 * t1214 + t1127 * t1203) * t1252 + (t1139 * t1099 + t1102 * t1226) * t1175 - t1142 * t1246 + t1123 * t1218) * t1193 - ((t1110 * t1215 + t1127 * t1204) * t1253 + (t1138 * t1098 + t1101 * t1227) * t1173 - t1141 * t1247 + t1122 * t1220) * t1194 - ((t1110 * t1216 + t1127 * t1205) * t1254 + (t1137 * t1097 + t1100 * t1228) * t1171 - t1140 * t1248 + t1121 * t1222) * t1195) * m(3); ((-t1127 * t1217 + t1206) * t1175 + (t1170 * t1203 + t1226) * t1169) * t1196 + t1142 * t1229 + ((-t1127 * t1219 + t1207) * t1173 + (t1168 * t1204 + t1227) * t1167) * t1197 + t1141 * t1230 + ((-t1127 * t1221 + t1208) * t1171 + (t1166 * t1205 + t1228) * t1165) * t1198 + t1140 * t1231 - m(4) * g(2) + (-((-t1110 * t1217 + t1127 * t1206) * t1252 + (t1099 * t1142 - t1102 * t1218) * t1175 + t1139 * t1246 + t1123 * t1226) * t1193 - ((-t1110 * t1219 + t1127 * t1207) * t1253 + (t1098 * t1141 - t1101 * t1220) * t1173 + t1138 * t1247 + t1122 * t1227) * t1194 - ((-t1110 * t1221 + t1127 * t1208) * t1254 + (t1097 * t1140 - t1100 * t1222) * t1171 + t1137 * t1248 + t1121 * t1228) * t1195) * m(3); -m(4) * g(3) + t1172 * t1225 + t1174 * t1224 + t1176 * t1223 + (-(t1170 * t1155 + (pkin(1) + t1189) * t1176) * t1238 - (t1168 * t1154 + (pkin(1) + t1190) * t1174) * t1239 - (t1166 * t1153 + (pkin(1) + t1191) * t1172) * t1240) * m(3);];
taugX  = t1;
