% Calculate Gravitation load for parallel robot
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:36
% EndTime: 2020-08-07 10:58:37
% DurationCPUTime: 0.74s
% Computational Cost: add. (479->128), mult. (1062->248), div. (128->13), fcn. (746->26), ass. (0->119)
t1171 = sin(qJ(3,1));
t1177 = cos(qJ(3,1));
t1237 = -rSges(3,1) * t1177 + rSges(3,2) * t1171;
t1169 = sin(qJ(3,2));
t1175 = cos(qJ(3,2));
t1236 = -rSges(3,1) * t1175 + rSges(3,2) * t1169;
t1167 = sin(qJ(3,3));
t1173 = cos(qJ(3,3));
t1235 = -rSges(3,1) * t1173 + rSges(3,2) * t1167;
t1159 = sin(qJ(3,4));
t1161 = cos(qJ(3,4));
t1234 = -rSges(3,1) * t1161 + rSges(3,2) * t1159;
t1233 = g(3) * rSges(3,3);
t1163 = legFrame(4,2);
t1136 = sin(t1163);
t1140 = cos(t1163);
t1132 = g(1) * t1140 - g(2) * t1136;
t1147 = 0.1e1 / t1161;
t1128 = g(1) * t1136 + g(2) * t1140;
t1160 = sin(qJ(2,4));
t1162 = cos(qJ(2,4));
t1204 = g(3) * t1162 + t1160 * t1128;
t1224 = ((rSges(3,1) * t1132 + t1204 * rSges(3,2)) * t1161 + t1159 * (t1204 * rSges(3,1) - rSges(3,2) * t1132)) * t1147;
t1164 = legFrame(3,2);
t1137 = sin(t1164);
t1141 = cos(t1164);
t1133 = g(1) * t1141 - g(2) * t1137;
t1153 = 0.1e1 / t1173;
t1129 = g(1) * t1137 + g(2) * t1141;
t1168 = sin(qJ(2,3));
t1174 = cos(qJ(2,3));
t1203 = g(3) * t1174 + t1168 * t1129;
t1223 = ((rSges(3,1) * t1133 + t1203 * rSges(3,2)) * t1173 + t1167 * (t1203 * rSges(3,1) - rSges(3,2) * t1133)) * t1153;
t1165 = legFrame(2,2);
t1138 = sin(t1165);
t1142 = cos(t1165);
t1134 = g(1) * t1142 - g(2) * t1138;
t1155 = 0.1e1 / t1175;
t1130 = g(1) * t1138 + g(2) * t1142;
t1170 = sin(qJ(2,2));
t1176 = cos(qJ(2,2));
t1202 = g(3) * t1176 + t1170 * t1130;
t1222 = ((rSges(3,1) * t1134 + t1202 * rSges(3,2)) * t1175 + t1169 * (t1202 * rSges(3,1) - rSges(3,2) * t1134)) * t1155;
t1166 = legFrame(1,2);
t1139 = sin(t1166);
t1143 = cos(t1166);
t1135 = g(1) * t1143 - g(2) * t1139;
t1157 = 0.1e1 / t1177;
t1131 = g(1) * t1139 + g(2) * t1143;
t1172 = sin(qJ(2,1));
t1178 = cos(qJ(2,1));
t1201 = g(3) * t1178 + t1172 * t1131;
t1221 = ((rSges(3,1) * t1135 + t1201 * rSges(3,2)) * t1177 + t1171 * (t1201 * rSges(3,1) - rSges(3,2) * t1135)) * t1157;
t1146 = 0.1e1 / t1160;
t1220 = t1146 * t1147;
t1219 = t1146 * t1162;
t1150 = 0.1e1 / t1168;
t1218 = t1150 * t1153;
t1217 = t1150 * t1174;
t1151 = 0.1e1 / t1170;
t1216 = t1151 * t1155;
t1215 = t1151 * t1176;
t1152 = 0.1e1 / t1172;
t1214 = t1152 * t1157;
t1213 = t1152 * t1178;
t1212 = t1160 * t1161;
t1211 = t1168 * t1173;
t1210 = t1170 * t1175;
t1209 = t1172 * t1177;
t1208 = t1128 * t1220;
t1207 = t1129 * t1218;
t1206 = t1130 * t1216;
t1205 = t1131 * t1214;
t1179 = rSges(2,2) * g(3);
t1180 = rSges(2,1) * g(3);
t1104 = ((-rSges(2,1) * t1128 + t1179) * t1162 + t1160 * (rSges(2,2) * t1128 + t1180)) * m(2) + ((t1234 * t1128 - t1233) * t1162 + t1160 * (-t1128 * rSges(3,3) - t1234 * g(3))) * m(3);
t1200 = 0.1e1 / t1161 ^ 2 * t1159 * t1104 * t1219;
t1105 = ((-rSges(2,1) * t1129 + t1179) * t1174 + t1168 * (rSges(2,2) * t1129 + t1180)) * m(2) + ((t1235 * t1129 - t1233) * t1174 + t1168 * (-t1129 * rSges(3,3) - t1235 * g(3))) * m(3);
t1199 = 0.1e1 / t1173 ^ 2 * t1167 * t1105 * t1217;
t1106 = ((-rSges(2,1) * t1130 + t1179) * t1176 + t1170 * (rSges(2,2) * t1130 + t1180)) * m(2) + ((t1236 * t1130 - t1233) * t1176 + t1170 * (-t1130 * rSges(3,3) - t1236 * g(3))) * m(3);
t1198 = 0.1e1 / t1175 ^ 2 * t1169 * t1106 * t1215;
t1107 = ((-rSges(2,1) * t1131 + t1179) * t1178 + t1172 * (rSges(2,2) * t1131 + t1180)) * m(2) + ((t1237 * t1131 - t1233) * t1178 + t1172 * (-t1131 * rSges(3,3) - t1237 * g(3))) * m(3);
t1197 = 0.1e1 / t1177 ^ 2 * t1171 * t1107 * t1213;
t1181 = xP(4);
t1144 = sin(t1181);
t1145 = cos(t1181);
t1184 = koppelP(4,2);
t1188 = koppelP(4,1);
t1120 = -t1144 * t1188 - t1145 * t1184;
t1124 = -t1144 * t1184 + t1145 * t1188;
t1196 = t1120 * t1140 - t1124 * t1136;
t1185 = koppelP(3,2);
t1189 = koppelP(3,1);
t1121 = -t1144 * t1189 - t1145 * t1185;
t1125 = -t1144 * t1185 + t1145 * t1189;
t1195 = t1121 * t1141 - t1125 * t1137;
t1186 = koppelP(2,2);
t1190 = koppelP(2,1);
t1122 = -t1144 * t1190 - t1145 * t1186;
t1126 = -t1144 * t1186 + t1145 * t1190;
t1194 = t1122 * t1142 - t1126 * t1138;
t1187 = koppelP(1,2);
t1191 = koppelP(1,1);
t1123 = -t1144 * t1191 - t1145 * t1187;
t1127 = -t1144 * t1187 + t1145 * t1191;
t1193 = t1123 * t1143 - t1127 * t1139;
t1192 = 0.1e1 / pkin(2);
t1183 = rSges(4,1);
t1182 = rSges(4,2);
t1149 = m(1) + m(2) + m(3);
t1119 = t1139 * t1171 + t1143 * t1209;
t1118 = t1138 * t1169 + t1142 * t1210;
t1117 = t1137 * t1167 + t1141 * t1211;
t1116 = t1139 * t1209 - t1143 * t1171;
t1115 = t1138 * t1210 - t1142 * t1169;
t1114 = t1137 * t1211 - t1141 * t1167;
t1113 = t1136 * t1159 + t1140 * t1212;
t1112 = t1136 * t1212 - t1140 * t1159;
t1 = [-m(4) * g(1) + (-t1112 * t1208 - t1114 * t1207 - t1115 * t1206 - t1116 * t1205) * t1149 + (t1140 * t1200 + t1141 * t1199 + t1142 * t1198 + t1143 * t1197 + (-t1140 * t1224 - t1141 * t1223 - t1142 * t1222 - t1143 * t1221) * m(3)) * t1192; -m(4) * g(2) + (-t1113 * t1208 - t1117 * t1207 - t1118 * t1206 - t1119 * t1205) * t1149 + (-t1136 * t1200 - t1137 * t1199 - t1138 * t1198 - t1139 * t1197 + (t1136 * t1224 + t1137 * t1223 + t1138 * t1222 + t1139 * t1221) * m(3)) * t1192; -m(4) * g(3) + (-t1104 * t1220 - t1105 * t1218 - t1106 * t1216 - t1107 * t1214) * t1192 + (-t1128 * t1219 - t1129 * t1217 - t1130 * t1215 - t1131 * t1213) * t1149; ((g(1) * t1182 - g(2) * t1183) * t1145 + t1144 * (g(1) * t1183 + g(2) * t1182)) * m(4) + (-(t1116 * t1123 + t1119 * t1127) * t1205 - (t1115 * t1122 + t1118 * t1126) * t1206 - (t1114 * t1121 + t1117 * t1125) * t1207 - (t1112 * t1120 + t1113 * t1124) * t1208) * t1149 + (t1193 * t1197 + t1194 * t1198 + t1195 * t1199 + t1196 * t1200 + (-t1193 * t1221 - t1194 * t1222 - t1195 * t1223 - t1196 * t1224) * m(3)) * t1192;];
taugX  = t1;
