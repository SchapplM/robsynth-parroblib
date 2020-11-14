% Calculate Gravitation load for parallel robot
% P3RPRRR9V1G3A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:56:26
% EndTime: 2020-08-06 18:56:28
% DurationCPUTime: 1.28s
% Computational Cost: add. (672->183), mult. (1095->290), div. (42->10), fcn. (672->26), ass. (0->130)
t1152 = m(2) + m(3);
t1179 = m(1) * rSges(1,1) + pkin(1) * t1152;
t1219 = 2 * pkin(3);
t1134 = cos(pkin(7));
t1218 = 0.2e1 * t1134 ^ 2;
t1217 = m(2) * rSges(2,1);
t1216 = (pkin(5) + qJ(2,1));
t1215 = (pkin(5) + qJ(2,2));
t1214 = (pkin(5) + qJ(2,3));
t1213 = m(3) / pkin(3);
t1144 = cos(qJ(3,3));
t1130 = t1144 ^ 2;
t1212 = t1130 * pkin(3);
t1146 = cos(qJ(3,2));
t1131 = t1146 ^ 2;
t1211 = t1131 * pkin(3);
t1148 = cos(qJ(3,1));
t1132 = t1148 ^ 2;
t1210 = t1132 * pkin(3);
t1209 = t1144 * pkin(2);
t1208 = t1146 * pkin(2);
t1207 = t1148 * pkin(2);
t1135 = legFrame(3,2);
t1111 = sin(t1135);
t1114 = cos(t1135);
t1074 = t1111 * g(1) + t1114 * g(2);
t1127 = pkin(7) + qJ(3,3);
t1104 = cos(t1127);
t1098 = 0.1e1 / t1104;
t1101 = sin(t1127);
t1077 = t1114 * g(1) - t1111 * g(2);
t1139 = sin(qJ(1,3));
t1145 = cos(qJ(1,3));
t1165 = -t1139 * g(3) + t1077 * t1145;
t1206 = ((-rSges(3,1) * t1074 + t1165 * rSges(3,2)) * t1104 + t1101 * (t1165 * rSges(3,1) + rSges(3,2) * t1074)) * t1098;
t1136 = legFrame(2,2);
t1112 = sin(t1136);
t1115 = cos(t1136);
t1075 = t1112 * g(1) + t1115 * g(2);
t1128 = pkin(7) + qJ(3,2);
t1105 = cos(t1128);
t1099 = 0.1e1 / t1105;
t1102 = sin(t1128);
t1078 = t1115 * g(1) - t1112 * g(2);
t1141 = sin(qJ(1,2));
t1147 = cos(qJ(1,2));
t1164 = -t1141 * g(3) + t1078 * t1147;
t1205 = ((-rSges(3,1) * t1075 + t1164 * rSges(3,2)) * t1105 + t1102 * (t1164 * rSges(3,1) + rSges(3,2) * t1075)) * t1099;
t1137 = legFrame(1,2);
t1113 = sin(t1137);
t1116 = cos(t1137);
t1076 = t1113 * g(1) + t1116 * g(2);
t1129 = pkin(7) + qJ(3,1);
t1106 = cos(t1129);
t1100 = 0.1e1 / t1106;
t1103 = sin(t1129);
t1079 = t1116 * g(1) - t1113 * g(2);
t1143 = sin(qJ(1,1));
t1149 = cos(qJ(1,1));
t1163 = -t1143 * g(3) + t1079 * t1149;
t1204 = ((-rSges(3,1) * t1076 + t1163 * rSges(3,2)) * t1106 + t1103 * (t1163 * rSges(3,1) + rSges(3,2) * t1076)) * t1100;
t1123 = pkin(6) + t1214;
t1108 = 0.1e1 / t1123;
t1203 = (-g(3) * t1145 - t1077 * t1139) * t1108;
t1124 = pkin(6) + t1215;
t1109 = 0.1e1 / t1124;
t1202 = (-g(3) * t1147 - t1078 * t1141) * t1109;
t1125 = pkin(6) + t1216;
t1110 = 0.1e1 / t1125;
t1201 = (-g(3) * t1149 - t1079 * t1143) * t1110;
t1138 = sin(qJ(3,3));
t1155 = pkin(2) / 0.2e1;
t1200 = (t1144 * pkin(3) + t1155) * t1138;
t1140 = sin(qJ(3,2));
t1199 = (t1146 * pkin(3) + t1155) * t1140;
t1142 = sin(qJ(3,1));
t1198 = (t1148 * pkin(3) + t1155) * t1142;
t1133 = sin(pkin(7));
t1150 = rSges(2,2) * m(2);
t1162 = -(-t1217 + (-rSges(3,1) * t1144 + rSges(3,2) * t1138 - pkin(2)) * m(3)) * t1134 - (t1150 + (rSges(3,1) * t1138 + rSges(3,2) * t1144) * m(3)) * t1133;
t1166 = t1179 * g(3);
t1151 = (m(1) * rSges(1,2));
t1169 = (rSges(3,3) + t1214) * m(3) - (-rSges(2,3) - qJ(2,3)) * m(2) - t1151;
t1197 = t1108 * (t1166 * t1145 + (t1139 * t1169 + t1162 * t1145) * g(3) + (-t1169 * t1145 + t1139 * (t1162 + t1179)) * t1077);
t1161 = -(-t1217 + (-rSges(3,1) * t1146 + rSges(3,2) * t1140 - pkin(2)) * m(3)) * t1134 - (t1150 + (rSges(3,1) * t1140 + rSges(3,2) * t1146) * m(3)) * t1133;
t1168 = (rSges(3,3) + t1215) * m(3) - (-rSges(2,3) - qJ(2,2)) * m(2) - t1151;
t1196 = t1109 * (t1166 * t1147 + (t1141 * t1168 + t1161 * t1147) * g(3) + (-t1168 * t1147 + t1141 * (t1161 + t1179)) * t1078);
t1160 = -(-t1217 + (-rSges(3,1) * t1148 + rSges(3,2) * t1142 - pkin(2)) * m(3)) * t1134 - (t1150 + (rSges(3,1) * t1142 + rSges(3,2) * t1148) * m(3)) * t1133;
t1167 = (rSges(3,3) + t1216) * m(3) - (-rSges(2,3) - qJ(2,1)) * m(2) - t1151;
t1195 = t1110 * (t1166 * t1149 + (t1143 * t1167 + t1160 * t1149) * g(3) + (-t1167 * t1149 + t1143 * (t1160 + t1179)) * t1079);
t1194 = t1111 * t1145;
t1193 = t1112 * t1147;
t1192 = t1113 * t1149;
t1191 = t1114 * t1145;
t1190 = t1115 * t1147;
t1189 = t1116 * t1149;
t1188 = t1133 * t1138;
t1187 = t1133 * t1140;
t1186 = t1133 * t1142;
t1107 = pkin(1) * t1133;
t1185 = t1144 * (-t1138 * pkin(3) + t1107);
t1184 = t1146 * (-t1140 * pkin(3) + t1107);
t1183 = t1148 * (-t1142 * pkin(3) + t1107);
t1182 = pkin(1) * t1145 + t1139 * t1123;
t1181 = pkin(1) * t1147 + t1141 * t1124;
t1180 = pkin(1) * t1149 + t1143 * t1125;
t1178 = t1098 * t1197;
t1177 = t1099 * t1196;
t1176 = t1100 * t1195;
t1175 = 0.1e1 / (t1134 * t1144 - t1188) * t1203;
t1174 = 0.1e1 / (t1134 * t1146 - t1187) * t1202;
t1173 = 0.1e1 / (t1134 * t1148 - t1186) * t1201;
t1172 = t1145 * t1188;
t1171 = t1147 * t1187;
t1170 = t1149 * t1186;
t1159 = pkin(2) * t1172 + (t1172 * t1219 - t1182) * t1144;
t1158 = pkin(2) * t1171 + (t1171 * t1219 - t1181) * t1146;
t1157 = pkin(2) * t1170 + (t1170 * t1219 - t1180) * t1148;
t1154 = -pkin(3) / 0.2e1;
t1096 = t1134 * pkin(2) + pkin(1);
t1082 = t1210 + t1207 / 0.2e1 + t1154;
t1081 = t1211 + t1208 / 0.2e1 + t1154;
t1080 = t1212 + t1209 / 0.2e1 + t1154;
t1058 = pkin(1) * t1142 + (-pkin(3) + t1207 + 0.2e1 * t1210) * t1133;
t1057 = pkin(1) * t1140 + (-pkin(3) + t1208 + 0.2e1 * t1211) * t1133;
t1056 = pkin(1) * t1138 + (-pkin(3) + t1209 + 0.2e1 * t1212) * t1133;
t1055 = t1180 * t1186 + (t1132 - 0.1e1) * t1149 * pkin(3);
t1054 = t1181 * t1187 + (t1131 - 0.1e1) * t1147 * pkin(3);
t1053 = t1182 * t1188 + (t1130 - 0.1e1) * t1145 * pkin(3);
t1 = [(t1113 * t1103 + t1106 * t1189) * t1176 + (t1112 * t1102 + t1105 * t1190) * t1177 + (t1111 * t1101 + t1104 * t1191) * t1178 - m(4) * g(1) + (((t1082 * t1189 + t1113 * t1198) * t1218 + (t1113 * t1058 - t1157 * t1116) * t1134 - t1055 * t1116 + t1113 * t1183) * t1173 + ((t1081 * t1190 + t1112 * t1199) * t1218 + (t1112 * t1057 - t1158 * t1115) * t1134 - t1054 * t1115 + t1112 * t1184) * t1174 + ((t1080 * t1191 + t1111 * t1200) * t1218 + (t1111 * t1056 - t1159 * t1114) * t1134 - t1053 * t1114 + t1111 * t1185) * t1175) * t1152 + (t1111 * t1206 + t1112 * t1205 + t1113 * t1204) * t1213; (t1116 * t1103 - t1106 * t1192) * t1176 + (t1115 * t1102 - t1105 * t1193) * t1177 + (t1114 * t1101 - t1104 * t1194) * t1178 - m(4) * g(2) + (((-t1082 * t1192 + t1116 * t1198) * t1218 + (t1116 * t1058 + t1157 * t1113) * t1134 + t1055 * t1113 + t1116 * t1183) * t1173 + ((-t1081 * t1193 + t1115 * t1199) * t1218 + (t1115 * t1057 + t1158 * t1112) * t1134 + t1054 * t1112 + t1115 * t1184) * t1174 + ((-t1080 * t1194 + t1114 * t1200) * t1218 + (t1114 * t1056 + t1159 * t1111) * t1134 + t1053 * t1111 + t1114 * t1185) * t1175) * t1152 + (t1114 * t1206 + t1115 * t1205 + t1116 * t1204) * t1213; -t1139 * t1197 - t1141 * t1196 - t1143 * t1195 - m(4) * g(3) + ((t1125 * t1149 + (-pkin(3) * t1106 - t1096) * t1143) * t1201 + (t1124 * t1147 + (-pkin(3) * t1105 - t1096) * t1141) * t1202 + (t1123 * t1145 + (-pkin(3) * t1104 - t1096) * t1139) * t1203) * t1152;];
taugX  = t1;
