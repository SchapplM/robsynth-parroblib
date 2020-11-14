% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR12V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x14]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR12V1G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:27
% EndTime: 2020-08-06 18:28:28
% DurationCPUTime: 0.78s
% Computational Cost: add. (405->123), mult. (750->232), div. (102->7), fcn. (738->18), ass. (0->108)
t1107 = legFrame(3,2);
t1094 = sin(t1107);
t1097 = cos(t1107);
t1082 = g(1) * t1097 - g(2) * t1094;
t1111 = sin(qJ(1,3));
t1117 = cos(qJ(1,3));
t1074 = -g(3) * t1111 + t1082 * t1117;
t1110 = sin(qJ(3,3));
t1088 = pkin(3) * t1110 + qJ(2,3);
t1085 = 0.1e1 / t1088;
t1170 = t1085 * t1111;
t1149 = t1074 * t1170;
t1108 = legFrame(2,2);
t1095 = sin(t1108);
t1098 = cos(t1108);
t1083 = g(1) * t1098 - g(2) * t1095;
t1113 = sin(qJ(1,2));
t1119 = cos(qJ(1,2));
t1076 = -g(3) * t1113 + t1083 * t1119;
t1112 = sin(qJ(3,2));
t1089 = pkin(3) * t1112 + qJ(2,2);
t1086 = 0.1e1 / t1089;
t1168 = t1086 * t1113;
t1147 = t1076 * t1168;
t1109 = legFrame(1,2);
t1096 = sin(t1109);
t1099 = cos(t1109);
t1084 = g(1) * t1099 - g(2) * t1096;
t1115 = sin(qJ(1,1));
t1121 = cos(qJ(1,1));
t1078 = -g(3) * t1115 + t1084 * t1121;
t1114 = sin(qJ(3,1));
t1090 = pkin(3) * t1114 + qJ(2,1);
t1087 = 0.1e1 / t1090;
t1166 = t1087 * t1115;
t1145 = t1078 * t1166;
t1077 = g(3) * t1121 + t1084 * t1115;
t1075 = g(3) * t1119 + t1083 * t1113;
t1073 = g(3) * t1117 + t1082 * t1111;
t1104 = 0.1e1 / t1110;
t1176 = t1073 * t1104;
t1105 = 0.1e1 / t1112;
t1175 = t1075 * t1105;
t1106 = 0.1e1 / t1114;
t1174 = t1077 * t1106;
t1169 = t1085 * t1117;
t1167 = t1086 * t1119;
t1165 = t1087 * t1121;
t1164 = t1094 * t1104;
t1116 = cos(qJ(3,3));
t1163 = t1094 * t1116;
t1162 = t1094 * t1117;
t1161 = t1095 * t1105;
t1118 = cos(qJ(3,2));
t1160 = t1095 * t1118;
t1159 = t1095 * t1119;
t1158 = t1096 * t1106;
t1120 = cos(qJ(3,1));
t1157 = t1096 * t1120;
t1156 = t1096 * t1121;
t1155 = t1097 * t1104;
t1154 = t1097 * t1116;
t1153 = t1098 * t1105;
t1152 = t1098 * t1118;
t1151 = t1099 * t1106;
t1150 = t1099 * t1120;
t1148 = t1074 * t1169;
t1146 = t1076 * t1167;
t1144 = t1078 * t1165;
t1143 = t1085 * t1162;
t1142 = t1097 * t1169;
t1141 = t1086 * t1159;
t1140 = t1098 * t1167;
t1139 = t1087 * t1156;
t1138 = t1099 * t1165;
t1091 = t1115 * qJ(2,1);
t1103 = pkin(1) + pkin(5) + pkin(6);
t1137 = t1103 * t1121 + t1091;
t1092 = qJ(2,3) * t1111;
t1136 = t1103 * t1117 + t1092;
t1093 = qJ(2,2) * t1113;
t1135 = t1103 * t1119 + t1093;
t1134 = t1110 * t1148;
t1133 = t1116 * t1148;
t1132 = t1112 * t1146;
t1131 = t1118 * t1146;
t1130 = t1114 * t1144;
t1129 = t1120 * t1144;
t1128 = t1073 * t1170 + t1075 * t1168 + t1077 * t1166;
t1127 = t1145 + t1147 + t1149;
t1126 = t1073 * t1143 + t1075 * t1141 + t1077 * t1139;
t1125 = t1073 * t1142 + t1075 * t1140 + t1077 * t1138;
t1124 = t1074 * t1143 + t1076 * t1141 + t1078 * t1139;
t1123 = t1074 * t1142 + t1076 * t1140 + t1078 * t1138;
t1122 = 0.1e1 / pkin(3);
t1081 = g(1) * t1096 + g(2) * t1099;
t1080 = g(1) * t1095 + g(2) * t1098;
t1079 = g(1) * t1094 + g(2) * t1097;
t1066 = t1084 * (pkin(1) * t1115 - qJ(2,1) * t1121) + g(3) * (pkin(1) * t1121 + t1091);
t1065 = t1083 * (pkin(1) * t1113 - qJ(2,2) * t1119) + g(3) * (pkin(1) * t1119 + t1093);
t1064 = t1082 * (pkin(1) * t1111 - qJ(2,3) * t1117) + g(3) * (pkin(1) * t1117 + t1092);
t1063 = t1077 * t1114 + t1081 * t1120;
t1062 = -t1077 * t1120 + t1081 * t1114;
t1061 = t1075 * t1112 + t1080 * t1118;
t1060 = -t1075 * t1118 + t1080 * t1112;
t1059 = t1073 * t1110 + t1079 * t1116;
t1058 = -t1073 * t1116 + t1079 * t1110;
t1 = [0, t1125, t1123, -t1125, -t1123, t1066 * t1138 - (t1137 * t1099 * t1114 + qJ(2,1) * t1157 + (t1114 * t1157 + (-t1120 ^ 2 + 0.1e1) * t1099 * t1115) * pkin(3)) * t1087 * t1174 + t1065 * t1140 - (t1135 * t1098 * t1112 + qJ(2,2) * t1160 + (t1112 * t1160 + (-t1118 ^ 2 + 0.1e1) * t1098 * t1113) * pkin(3)) * t1086 * t1175 + t1064 * t1142 - (t1136 * t1097 * t1110 + qJ(2,3) * t1163 + (t1110 * t1163 + (-t1116 ^ 2 + 0.1e1) * t1097 * t1111) * pkin(3)) * t1085 * t1176, 0, 0, 0, 0, 0, -t1097 * t1134 - t1098 * t1132 - t1099 * t1130 + (-t1058 * t1164 - t1060 * t1161 - t1062 * t1158) * t1122, -t1097 * t1133 - t1098 * t1131 - t1099 * t1129 + (-t1059 * t1164 - t1061 * t1161 - t1063 * t1158) * t1122, -g(1); 0, -t1126, -t1124, t1126, t1124, (-t1066 * t1156 - ((pkin(3) * t1150 - t1137 * t1096) * t1114 + t1115 * pkin(3) * (t1120 - 0.1e1) * (t1120 + 0.1e1) * t1096 + qJ(2,1) * t1150) * t1174) * t1087 + (-t1065 * t1159 - ((pkin(3) * t1152 - t1135 * t1095) * t1112 + t1113 * pkin(3) * (t1118 - 0.1e1) * (t1118 + 0.1e1) * t1095 + qJ(2,2) * t1152) * t1175) * t1086 + (-t1064 * t1162 - ((pkin(3) * t1154 - t1136 * t1094) * t1110 + t1111 * pkin(3) * (t1116 - 0.1e1) * (t1116 + 0.1e1) * t1094 + qJ(2,3) * t1154) * t1176) * t1085, 0, 0, 0, 0, 0, t1094 * t1134 + t1095 * t1132 + t1096 * t1130 + (-t1058 * t1155 - t1060 * t1153 - t1062 * t1151) * t1122, t1094 * t1133 + t1095 * t1131 + t1096 * t1129 + (-t1059 * t1155 - t1061 * t1153 - t1063 * t1151) * t1122, -g(2); 0, -t1128, -t1127, t1128, t1127, (-t1115 * t1066 - (t1090 * t1121 - t1103 * t1115) * t1077) * t1087 + (-t1113 * t1065 - (t1089 * t1119 - t1103 * t1113) * t1075) * t1086 + (-t1111 * t1064 - (t1088 * t1117 - t1103 * t1111) * t1073) * t1085, 0, 0, 0, 0, 0, t1110 * t1149 + t1112 * t1147 + t1114 * t1145, t1116 * t1149 + t1118 * t1147 + t1120 * t1145, -g(3);];
tau_reg  = t1;
