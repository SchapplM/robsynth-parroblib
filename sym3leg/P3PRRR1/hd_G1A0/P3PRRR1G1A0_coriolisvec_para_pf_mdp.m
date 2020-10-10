% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:58
% EndTime: 2020-03-09 21:15:00
% DurationCPUTime: 1.63s
% Computational Cost: add. (14584->151), mult. (10130->307), div. (1608->9), fcn. (10128->24), ass. (0->156)
t1115 = pkin(7) + qJ(2,3);
t1106 = qJ(3,3) + t1115;
t1094 = sin(t1106);
t1097 = cos(t1106);
t1100 = sin(t1115);
t1103 = cos(t1115);
t1073 = t1094 * t1103 - t1100 * t1097;
t1212 = 0.1e1 / t1073;
t1116 = pkin(7) + qJ(2,2);
t1107 = qJ(3,2) + t1116;
t1095 = sin(t1107);
t1098 = cos(t1107);
t1101 = sin(t1116);
t1104 = cos(t1116);
t1074 = t1095 * t1104 - t1101 * t1098;
t1211 = 0.1e1 / t1074;
t1117 = pkin(7) + qJ(2,1);
t1108 = qJ(3,1) + t1117;
t1096 = sin(t1108);
t1099 = cos(t1108);
t1102 = sin(t1117);
t1105 = cos(t1117);
t1075 = t1096 * t1105 - t1102 * t1099;
t1210 = 0.1e1 / t1075;
t1065 = 0.1e1 / t1073 ^ 2;
t1068 = 0.1e1 / t1074 ^ 2;
t1071 = 0.1e1 / t1075 ^ 2;
t1131 = 0.1e1 / pkin(2);
t1130 = 0.1e1 / pkin(3);
t1118 = legFrame(3,3);
t1109 = sin(t1118);
t1112 = cos(t1118);
t1127 = xDP(2);
t1128 = xDP(1);
t1085 = t1109 * t1127 + t1112 * t1128;
t1086 = t1109 * t1128 - t1112 * t1127;
t1052 = -t1085 * t1097 + t1086 * t1094;
t1200 = (-(t1085 * t1103 - t1100 * t1086) * pkin(2) + t1052 * pkin(3)) * t1212;
t1179 = t1130 * t1200;
t1194 = t1052 * t1212;
t1041 = (t1179 - t1194) * t1131;
t1144 = t1094 * t1100 + t1097 * t1103;
t1031 = pkin(3) * t1041 - t1144 * t1194;
t1180 = t1212 * t1200;
t1155 = t1041 * t1180;
t1132 = 0.1e1 / pkin(2) ^ 2;
t1182 = t1130 * t1132;
t1034 = (t1144 * pkin(2) + pkin(3)) * t1155 * t1182;
t1193 = t1052 * t1065;
t1129 = pkin(3) ^ 2;
t1181 = 0.2e1 * pkin(3) * t1131;
t1203 = (-t1041 * t1129 + (t1194 - t1144 * (-t1194 + t1179 / 0.2e1) * t1181) * pkin(2)) * t1130;
t1022 = -t1034 + (t1155 + (-t1031 - t1203) * t1193) * t1132;
t1209 = t1022 * t1212;
t1119 = legFrame(2,3);
t1110 = sin(t1119);
t1113 = cos(t1119);
t1087 = t1110 * t1127 + t1113 * t1128;
t1088 = t1110 * t1128 - t1113 * t1127;
t1053 = -t1087 * t1098 + t1088 * t1095;
t1199 = (-(t1087 * t1104 - t1101 * t1088) * pkin(2) + t1053 * pkin(3)) * t1211;
t1177 = t1130 * t1199;
t1192 = t1053 * t1211;
t1043 = (t1177 - t1192) * t1131;
t1142 = t1095 * t1101 + t1098 * t1104;
t1032 = pkin(3) * t1043 - t1142 * t1192;
t1178 = t1211 * t1199;
t1153 = t1043 * t1178;
t1035 = (t1142 * pkin(2) + pkin(3)) * t1153 * t1182;
t1191 = t1053 * t1068;
t1202 = (-t1043 * t1129 + (t1192 - t1142 * (-t1192 + t1177 / 0.2e1) * t1181) * pkin(2)) * t1130;
t1023 = -t1035 + (t1153 + (-t1032 - t1202) * t1191) * t1132;
t1208 = t1023 * t1211;
t1120 = legFrame(1,3);
t1111 = sin(t1120);
t1114 = cos(t1120);
t1089 = t1111 * t1127 + t1114 * t1128;
t1090 = t1111 * t1128 - t1114 * t1127;
t1054 = -t1089 * t1099 + t1090 * t1096;
t1198 = (-(t1089 * t1105 - t1102 * t1090) * pkin(2) + t1054 * pkin(3)) * t1210;
t1175 = t1130 * t1198;
t1190 = t1054 * t1210;
t1045 = (t1175 - t1190) * t1131;
t1140 = t1096 * t1102 + t1099 * t1105;
t1033 = t1045 * pkin(3) - t1140 * t1190;
t1176 = t1210 * t1198;
t1151 = t1045 * t1176;
t1036 = (t1140 * pkin(2) + pkin(3)) * t1151 * t1182;
t1189 = t1054 * t1071;
t1201 = (-t1045 * t1129 + (t1190 - t1140 * (-t1190 + t1175 / 0.2e1) * t1181) * pkin(2)) * t1130;
t1024 = -t1036 + (t1151 + (-t1033 - t1201) * t1189) * t1132;
t1207 = t1024 * t1210;
t1025 = (-t1031 * t1193 + t1155) * t1132;
t1206 = t1025 * t1212;
t1026 = (-t1032 * t1191 + t1153) * t1132;
t1205 = t1026 * t1211;
t1027 = (-t1033 * t1189 + t1151) * t1132;
t1204 = t1027 * t1210;
t1197 = t1052 ^ 2 * t1212 * t1065;
t1196 = t1053 ^ 2 * t1211 * t1068;
t1195 = t1054 ^ 2 * t1210 * t1071;
t1143 = t1112 * t1094 + t1109 * t1097;
t1188 = t1212 * t1143;
t1080 = t1094 * t1109 - t1112 * t1097;
t1187 = t1212 * t1080;
t1141 = t1113 * t1095 + t1110 * t1098;
t1186 = t1211 * t1141;
t1082 = t1095 * t1110 - t1113 * t1098;
t1185 = t1211 * t1082;
t1139 = t1114 * t1096 + t1111 * t1099;
t1184 = t1210 * t1139;
t1084 = t1096 * t1111 - t1114 * t1099;
t1183 = t1210 * t1084;
t1121 = sin(qJ(3,3));
t1174 = t1121 * t1197;
t1124 = cos(qJ(3,3));
t1173 = t1124 * t1197;
t1122 = sin(qJ(3,2));
t1172 = t1122 * t1196;
t1125 = cos(qJ(3,2));
t1171 = t1125 * t1196;
t1123 = sin(qJ(3,1));
t1170 = t1123 * t1195;
t1126 = cos(qJ(3,1));
t1169 = t1126 * t1195;
t1168 = t1121 * t1206;
t1167 = t1124 * t1206;
t1166 = t1122 * t1205;
t1165 = t1125 * t1205;
t1164 = t1123 * t1204;
t1163 = t1126 * t1204;
t1019 = -t1034 + (0.2e1 * t1155 + (-0.2e1 * t1031 - t1203) * t1193) * t1132;
t1162 = t1019 * t1188;
t1161 = t1019 * t1187;
t1020 = -t1035 + (0.2e1 * t1153 + (-0.2e1 * t1032 - t1202) * t1191) * t1132;
t1160 = t1020 * t1186;
t1159 = t1020 * t1185;
t1021 = -t1036 + (0.2e1 * t1151 + (-0.2e1 * t1033 - t1201) * t1189) * t1132;
t1158 = t1021 * t1184;
t1157 = t1021 * t1183;
t1156 = (t1179 - 0.2e1 * t1194) * t1131 * t1180;
t1154 = (t1177 - 0.2e1 * t1192) * t1131 * t1178;
t1152 = (t1175 - 0.2e1 * t1190) * t1131 * t1176;
t1150 = t1121 * t1156;
t1149 = t1124 * t1156;
t1148 = t1122 * t1154;
t1147 = t1125 * t1154;
t1146 = t1123 * t1152;
t1145 = t1126 * t1152;
t1060 = -pkin(2) * (t1102 * t1111 - t1114 * t1105) - t1084 * pkin(3);
t1059 = -pkin(2) * (t1101 * t1110 - t1113 * t1104) - t1082 * pkin(3);
t1058 = -pkin(2) * (t1100 * t1109 - t1112 * t1103) - t1080 * pkin(3);
t1057 = pkin(2) * (t1102 * t1114 + t1111 * t1105) + t1139 * pkin(3);
t1056 = pkin(2) * (t1101 * t1113 + t1110 * t1104) + t1141 * pkin(3);
t1055 = pkin(2) * (t1100 * t1112 + t1109 * t1103) + t1143 * pkin(3);
t1 = [(-t1124 * t1161 - t1125 * t1159 - t1126 * t1157) * MDP(6) + (t1121 * t1161 + t1122 * t1159 + t1123 * t1157) * MDP(7) + ((-t1025 * t1187 - t1026 * t1185 - t1027 * t1183) * MDP(2) + (-t1022 * t1187 - t1023 * t1185 - t1024 * t1183) * MDP(5)) * t1131 + ((-t1058 * t1167 - t1059 * t1165 - t1060 * t1163) * MDP(6) + (t1058 * t1168 + t1059 * t1166 + t1060 * t1164) * MDP(7) + ((-t1058 * t1174 - t1059 * t1172 - t1060 * t1170) * MDP(6) + (-t1058 * t1173 - t1059 * t1171 - t1060 * t1169) * MDP(7)) * t1132 + ((-t1058 * t1209 - t1059 * t1208 - t1060 * t1207) * MDP(5) + (t1080 * t1150 + t1082 * t1148 + t1084 * t1146) * MDP(6) + (t1080 * t1149 + t1082 * t1147 + t1084 * t1145) * MDP(7)) * t1131) * t1130; (t1124 * t1162 + t1125 * t1160 + t1126 * t1158) * MDP(6) + (-t1121 * t1162 - t1122 * t1160 - t1123 * t1158) * MDP(7) + ((t1025 * t1188 + t1026 * t1186 + t1027 * t1184) * MDP(2) + (t1022 * t1188 + t1023 * t1186 + t1024 * t1184) * MDP(5)) * t1131 + ((-t1055 * t1167 - t1056 * t1165 - t1057 * t1163) * MDP(6) + (t1055 * t1168 + t1056 * t1166 + t1057 * t1164) * MDP(7) + ((-t1055 * t1174 - t1056 * t1172 - t1057 * t1170) * MDP(6) + (-t1055 * t1173 - t1056 * t1171 - t1057 * t1169) * MDP(7)) * t1132 + ((-t1055 * t1209 - t1056 * t1208 - t1057 * t1207) * MDP(5) + (-t1139 * t1146 - t1141 * t1148 - t1143 * t1150) * MDP(6) + (-t1139 * t1145 - t1141 * t1147 - t1143 * t1149) * MDP(7)) * t1131) * t1130; 0;];
taucX  = t1;
