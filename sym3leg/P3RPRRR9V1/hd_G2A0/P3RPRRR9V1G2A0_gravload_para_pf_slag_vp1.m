% Calculate Gravitation load for parallel robot
% P3RPRRR9V1G2A0
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
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:28
% EndTime: 2020-08-06 18:51:29
% DurationCPUTime: 1.30s
% Computational Cost: add. (672->183), mult. (1095->290), div. (42->10), fcn. (672->26), ass. (0->130)
t1140 = m(2) + m(3);
t1209 = m(1) * rSges(1,1) + pkin(1) * t1140;
t1207 = 2 * pkin(3);
t1122 = cos(pkin(7));
t1206 = 0.2e1 * t1122 ^ 2;
t1205 = m(2) * rSges(2,1);
t1204 = (pkin(5) + qJ(2,1));
t1203 = (pkin(5) + qJ(2,2));
t1202 = (pkin(5) + qJ(2,3));
t1201 = m(3) / pkin(3);
t1132 = cos(qJ(3,3));
t1118 = t1132 ^ 2;
t1200 = t1118 * pkin(3);
t1134 = cos(qJ(3,2));
t1119 = t1134 ^ 2;
t1199 = t1119 * pkin(3);
t1136 = cos(qJ(3,1));
t1120 = t1136 ^ 2;
t1198 = t1120 * pkin(3);
t1197 = t1132 * pkin(2);
t1196 = t1134 * pkin(2);
t1195 = t1136 * pkin(2);
t1123 = legFrame(3,2);
t1099 = sin(t1123);
t1102 = cos(t1123);
t1065 = t1099 * g(1) + t1102 * g(2);
t1115 = pkin(7) + qJ(3,3);
t1092 = cos(t1115);
t1086 = 0.1e1 / t1092;
t1089 = sin(t1115);
t1068 = t1102 * g(1) - t1099 * g(2);
t1127 = sin(qJ(1,3));
t1133 = cos(qJ(1,3));
t1153 = g(3) * t1133 + t1068 * t1127;
t1194 = ((-rSges(3,1) * t1065 + t1153 * rSges(3,2)) * t1092 + t1089 * (t1153 * rSges(3,1) + rSges(3,2) * t1065)) * t1086;
t1124 = legFrame(2,2);
t1100 = sin(t1124);
t1103 = cos(t1124);
t1066 = t1100 * g(1) + t1103 * g(2);
t1116 = pkin(7) + qJ(3,2);
t1093 = cos(t1116);
t1087 = 0.1e1 / t1093;
t1090 = sin(t1116);
t1069 = t1103 * g(1) - t1100 * g(2);
t1129 = sin(qJ(1,2));
t1135 = cos(qJ(1,2));
t1152 = g(3) * t1135 + t1069 * t1129;
t1193 = ((-rSges(3,1) * t1066 + t1152 * rSges(3,2)) * t1093 + t1090 * (t1152 * rSges(3,1) + rSges(3,2) * t1066)) * t1087;
t1125 = legFrame(1,2);
t1101 = sin(t1125);
t1104 = cos(t1125);
t1067 = t1101 * g(1) + t1104 * g(2);
t1117 = pkin(7) + qJ(3,1);
t1094 = cos(t1117);
t1088 = 0.1e1 / t1094;
t1091 = sin(t1117);
t1070 = t1104 * g(1) - t1101 * g(2);
t1131 = sin(qJ(1,1));
t1137 = cos(qJ(1,1));
t1151 = g(3) * t1137 + t1070 * t1131;
t1192 = ((-rSges(3,1) * t1067 + t1151 * rSges(3,2)) * t1094 + t1091 * (t1151 * rSges(3,1) + rSges(3,2) * t1067)) * t1088;
t1111 = pkin(6) + t1202;
t1096 = 0.1e1 / t1111;
t1191 = (-g(3) * t1127 + t1068 * t1133) * t1096;
t1112 = pkin(6) + t1203;
t1097 = 0.1e1 / t1112;
t1190 = (-g(3) * t1129 + t1069 * t1135) * t1097;
t1113 = pkin(6) + t1204;
t1098 = 0.1e1 / t1113;
t1189 = (-g(3) * t1131 + t1070 * t1137) * t1098;
t1126 = sin(qJ(3,3));
t1143 = pkin(2) / 0.2e1;
t1188 = (t1132 * pkin(3) + t1143) * t1126;
t1128 = sin(qJ(3,2));
t1187 = (t1134 * pkin(3) + t1143) * t1128;
t1130 = sin(qJ(3,1));
t1186 = (t1136 * pkin(3) + t1143) * t1130;
t1121 = sin(pkin(7));
t1138 = rSges(2,2) * m(2);
t1150 = -(-t1205 + (-rSges(3,1) * t1132 + rSges(3,2) * t1126 - pkin(2)) * m(3)) * t1122 - (t1138 + (rSges(3,1) * t1126 + rSges(3,2) * t1132) * m(3)) * t1121;
t1157 = t1209 * g(3);
t1139 = (m(1) * rSges(1,2));
t1160 = -(rSges(3,3) + t1202) * m(3) + (-rSges(2,3) - qJ(2,3)) * m(2) + t1139;
t1185 = t1096 * (t1157 * t1127 + (t1150 * t1127 + t1160 * t1133) * g(3) + ((-t1150 - t1209) * t1133 + t1160 * t1127) * t1068);
t1149 = -(-t1205 + (-rSges(3,1) * t1134 + rSges(3,2) * t1128 - pkin(2)) * m(3)) * t1122 - (t1138 + (rSges(3,1) * t1128 + rSges(3,2) * t1134) * m(3)) * t1121;
t1159 = -(rSges(3,3) + t1203) * m(3) + (-rSges(2,3) - qJ(2,2)) * m(2) + t1139;
t1184 = t1097 * (t1157 * t1129 + (t1149 * t1129 + t1159 * t1135) * g(3) + ((-t1149 - t1209) * t1135 + t1159 * t1129) * t1069);
t1148 = -(-t1205 + (-rSges(3,1) * t1136 + rSges(3,2) * t1130 - pkin(2)) * m(3)) * t1122 - (t1138 + (rSges(3,1) * t1130 + rSges(3,2) * t1136) * m(3)) * t1121;
t1158 = -(rSges(3,3) + t1204) * m(3) + (-rSges(2,3) - qJ(2,1)) * m(2) + t1139;
t1183 = t1098 * (t1157 * t1131 + (t1148 * t1131 + t1158 * t1137) * g(3) + ((-t1148 - t1209) * t1137 + t1158 * t1131) * t1070);
t1182 = t1099 * t1127;
t1181 = t1100 * t1129;
t1180 = t1101 * t1131;
t1179 = t1102 * t1127;
t1178 = t1103 * t1129;
t1177 = t1104 * t1131;
t1176 = t1121 * t1126;
t1175 = t1121 * t1128;
t1174 = t1121 * t1130;
t1095 = pkin(1) * t1121;
t1173 = t1132 * (-t1126 * pkin(3) + t1095);
t1172 = t1134 * (-t1128 * pkin(3) + t1095);
t1171 = t1136 * (-t1130 * pkin(3) + t1095);
t1169 = t1086 * t1185;
t1168 = t1087 * t1184;
t1167 = t1088 * t1183;
t1166 = 0.1e1 / (t1122 * t1132 - t1176) * t1191;
t1165 = 0.1e1 / (t1122 * t1134 - t1175) * t1190;
t1164 = 0.1e1 / (t1122 * t1136 - t1174) * t1189;
t1163 = t1127 * t1176;
t1162 = t1129 * t1175;
t1161 = t1131 * t1174;
t1156 = pkin(1) * t1127 - t1133 * t1111;
t1155 = pkin(1) * t1129 - t1135 * t1112;
t1154 = pkin(1) * t1131 - t1137 * t1113;
t1147 = pkin(2) * t1163 + (t1163 * t1207 - t1156) * t1132;
t1146 = pkin(2) * t1162 + (t1162 * t1207 - t1155) * t1134;
t1145 = pkin(2) * t1161 + (t1161 * t1207 - t1154) * t1136;
t1142 = -pkin(3) / 0.2e1;
t1084 = t1122 * pkin(2) + pkin(1);
t1073 = t1198 + t1195 / 0.2e1 + t1142;
t1072 = t1199 + t1196 / 0.2e1 + t1142;
t1071 = t1200 + t1197 / 0.2e1 + t1142;
t1049 = pkin(1) * t1130 + (-pkin(3) + t1195 + 0.2e1 * t1198) * t1121;
t1048 = pkin(1) * t1128 + (-pkin(3) + t1196 + 0.2e1 * t1199) * t1121;
t1047 = pkin(1) * t1126 + (-pkin(3) + t1197 + 0.2e1 * t1200) * t1121;
t1046 = t1154 * t1174 + (t1120 - 0.1e1) * t1131 * pkin(3);
t1045 = t1155 * t1175 + (t1119 - 0.1e1) * t1129 * pkin(3);
t1044 = t1156 * t1176 + (t1118 - 0.1e1) * t1127 * pkin(3);
t1 = [(t1101 * t1091 + t1094 * t1177) * t1167 + (t1100 * t1090 + t1093 * t1178) * t1168 + (t1099 * t1089 + t1092 * t1179) * t1169 - m(4) * g(1) + (((t1073 * t1177 + t1101 * t1186) * t1206 + (t1101 * t1049 - t1104 * t1145) * t1122 - t1046 * t1104 + t1101 * t1171) * t1164 + ((t1072 * t1178 + t1100 * t1187) * t1206 + (t1100 * t1048 - t1103 * t1146) * t1122 - t1045 * t1103 + t1100 * t1172) * t1165 + ((t1071 * t1179 + t1099 * t1188) * t1206 + (t1099 * t1047 - t1102 * t1147) * t1122 - t1044 * t1102 + t1099 * t1173) * t1166) * t1140 + (t1099 * t1194 + t1100 * t1193 + t1101 * t1192) * t1201; (t1104 * t1091 - t1094 * t1180) * t1167 + (t1103 * t1090 - t1093 * t1181) * t1168 + (t1102 * t1089 - t1092 * t1182) * t1169 - m(4) * g(2) + (((-t1073 * t1180 + t1104 * t1186) * t1206 + (t1104 * t1049 + t1101 * t1145) * t1122 + t1046 * t1101 + t1104 * t1171) * t1164 + ((-t1072 * t1181 + t1103 * t1187) * t1206 + (t1103 * t1048 + t1100 * t1146) * t1122 + t1045 * t1100 + t1103 * t1172) * t1165 + ((-t1071 * t1182 + t1102 * t1188) * t1206 + (t1102 * t1047 + t1099 * t1147) * t1122 + t1044 * t1099 + t1102 * t1173) * t1166) * t1140 + (t1102 * t1194 + t1103 * t1193 + t1104 * t1192) * t1201; t1133 * t1185 + t1135 * t1184 + t1137 * t1183 - m(4) * g(3) + ((t1131 * t1113 + (pkin(3) * t1094 + t1084) * t1137) * t1189 + (t1129 * t1112 + (pkin(3) * t1093 + t1084) * t1135) * t1190 + (t1127 * t1111 + (pkin(3) * t1092 + t1084) * t1133) * t1191) * t1140;];
taugX  = t1;
