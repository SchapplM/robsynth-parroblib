% Calculate Gravitation load for parallel robot
% P3RRPRR12V2G3A0
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:27:34
% EndTime: 2020-08-06 19:27:35
% DurationCPUTime: 1.06s
% Computational Cost: add. (870->178), mult. (1047->312), div. (45->6), fcn. (648->18), ass. (0->131)
t1158 = mrSges(3,3) - mrSges(2,2);
t1052 = qJ(3,1) * m(3) + t1158;
t1053 = m(3) * pkin(2) + mrSges(2,1) + mrSges(3,1);
t1082 = sin(qJ(2,1));
t1088 = cos(qJ(2,1));
t1163 = t1082 * t1052 + t1053 * t1088;
t1051 = m(3) * qJ(3,2) + t1158;
t1080 = sin(qJ(2,2));
t1086 = cos(qJ(2,2));
t1162 = t1080 * t1051 + t1053 * t1086;
t1050 = m(3) * qJ(3,3) + t1158;
t1078 = sin(qJ(2,3));
t1084 = cos(qJ(2,3));
t1161 = t1078 * t1050 + t1053 * t1084;
t1091 = (pkin(2) + pkin(3));
t1160 = -2 * t1091;
t1159 = m(2) + m(3);
t1075 = legFrame(3,2);
t1057 = sin(t1075);
t1157 = t1057 * qJ(3,3);
t1076 = legFrame(2,2);
t1058 = sin(t1076);
t1156 = t1058 * qJ(3,2);
t1077 = legFrame(1,2);
t1059 = sin(t1077);
t1155 = t1059 * qJ(3,1);
t1060 = cos(t1075);
t1154 = t1060 * qJ(3,3);
t1061 = cos(t1076);
t1153 = t1061 * qJ(3,2);
t1062 = cos(t1077);
t1152 = t1062 * qJ(3,1);
t1054 = t1078 * qJ(3,3);
t1055 = t1080 * qJ(3,2);
t1056 = t1082 * qJ(3,1);
t1034 = t1057 * g(1) + t1060 * g(2);
t1037 = t1060 * g(1) - t1057 * g(2);
t1079 = sin(qJ(1,3));
t1085 = cos(qJ(1,3));
t1097 = g(3) * t1079 - t1037 * t1085;
t1013 = (-t1034 * t1053 + t1050 * t1097) * t1084 + t1078 * (-t1034 * t1050 - t1097 * t1053);
t1092 = 0.1e1 / qJ(3,3);
t1151 = t1013 * t1092;
t1035 = t1058 * g(1) + t1061 * g(2);
t1038 = t1061 * g(1) - t1058 * g(2);
t1081 = sin(qJ(1,2));
t1087 = cos(qJ(1,2));
t1096 = g(3) * t1081 - t1038 * t1087;
t1014 = (-t1035 * t1053 + t1051 * t1096) * t1086 + t1080 * (-t1035 * t1051 - t1096 * t1053);
t1093 = 0.1e1 / qJ(3,2);
t1150 = t1014 * t1093;
t1036 = t1059 * g(1) + t1062 * g(2);
t1039 = t1062 * g(1) - t1059 * g(2);
t1083 = sin(qJ(1,1));
t1089 = cos(qJ(1,1));
t1095 = g(3) * t1083 - t1039 * t1089;
t1015 = (-t1036 * t1053 + t1052 * t1095) * t1088 + t1082 * (-t1036 * t1052 - t1095 * t1053);
t1094 = 0.1e1 / qJ(3,1);
t1149 = t1015 * t1094;
t1107 = t1085 * t1054;
t1090 = pkin(5) - pkin(6);
t1115 = pkin(1) * t1085 + t1079 * t1090;
t1148 = (0.2e1 * t1107 + t1115) * t1091;
t1108 = t1087 * t1055;
t1114 = pkin(1) * t1087 + t1081 * t1090;
t1147 = (0.2e1 * t1108 + t1114) * t1091;
t1109 = t1089 * t1056;
t1113 = pkin(1) * t1089 + t1083 * t1090;
t1146 = (0.2e1 * t1109 + t1113) * t1091;
t1118 = t1091 * t1084;
t1100 = t1054 + pkin(1) + t1118;
t1031 = 0.1e1 / t1100;
t1145 = t1031 * t1092;
t1117 = t1091 * t1086;
t1099 = t1055 + pkin(1) + t1117;
t1032 = 0.1e1 / t1099;
t1144 = t1032 * t1093;
t1116 = t1091 * t1088;
t1098 = t1056 + pkin(1) + t1116;
t1033 = 0.1e1 / t1098;
t1143 = t1033 * t1094;
t1139 = (qJ(3,3) + t1091) * (-qJ(3,3) + t1091);
t1138 = (qJ(3,2) + t1091) * (-qJ(3,2) + t1091);
t1137 = (qJ(3,1) + t1091) * (-qJ(3,1) + t1091);
t1135 = t1078 * t1091;
t1041 = t1159 * pkin(5) - mrSges(1,2) + mrSges(3,2) + mrSges(2,3);
t1040 = g(3) * t1041;
t1043 = t1159 * pkin(1) + mrSges(1,1);
t1042 = g(3) * t1043;
t1010 = (t1161 * g(3) + t1042) * t1085 + t1040 * t1079 + (-t1041 * t1085 + (t1043 + t1161) * t1079) * t1037;
t1134 = t1079 * t1010;
t1132 = t1080 * t1091;
t1011 = (t1162 * g(3) + t1042) * t1087 + t1040 * t1081 + (-t1041 * t1087 + (t1043 + t1162) * t1081) * t1038;
t1131 = t1081 * t1011;
t1129 = t1082 * t1091;
t1012 = (t1163 * g(3) + t1042) * t1089 + t1040 * t1083 + (-t1041 * t1089 + (t1043 + t1163) * t1083) * t1039;
t1128 = t1083 * t1012;
t1127 = t1085 * t1091;
t1126 = t1087 * t1091;
t1125 = t1089 * t1091;
t1124 = t1090 * t1085;
t1123 = t1090 * t1087;
t1122 = t1090 * t1089;
t1044 = pkin(1) * t1078 + qJ(3,3);
t1121 = t1091 * t1044;
t1045 = pkin(1) * t1080 + qJ(3,2);
t1120 = t1091 * t1045;
t1046 = pkin(1) * t1082 + qJ(3,1);
t1119 = t1091 * t1046;
t1112 = qJ(3,1) * t1160;
t1111 = qJ(3,2) * t1160;
t1110 = qJ(3,3) * t1160;
t1106 = (-t1084 * t1034 - t1078 * t1097) * t1145;
t1105 = (-t1086 * t1035 - t1080 * t1096) * t1144;
t1104 = (-t1088 * t1036 - t1082 * t1095) * t1143;
t1103 = t1085 * t1139;
t1102 = t1087 * t1138;
t1101 = t1089 * t1137;
t1074 = t1088 ^ 2;
t1073 = t1086 ^ 2;
t1072 = t1084 ^ 2;
t1030 = pkin(1) * qJ(3,1) - t1082 * t1137;
t1029 = pkin(1) * qJ(3,2) - t1080 * t1138;
t1028 = pkin(1) * qJ(3,3) - t1078 * t1139;
t1027 = t1109 + t1113;
t1025 = t1108 + t1114;
t1023 = t1107 + t1115;
t1021 = qJ(3,1) * t1089 + t1082 * t1113;
t1020 = qJ(3,2) * t1087 + t1080 * t1114;
t1019 = qJ(3,3) * t1085 + t1078 * t1115;
t1 = [-g(1) * m(4) + (-t1062 * t1128 + ((t1062 * t1125 - t1155) * t1074 + (t1027 * t1062 + t1059 * t1129) * t1088 + t1059 * t1046) * t1149) * t1033 + (-t1061 * t1131 + ((t1061 * t1126 - t1156) * t1073 + (t1025 * t1061 + t1058 * t1132) * t1086 + t1058 * t1045) * t1150) * t1032 + (-t1060 * t1134 + ((t1060 * t1127 - t1157) * t1072 + (t1023 * t1060 + t1057 * t1135) * t1084 + t1057 * t1044) * t1151) * t1031 + (-((t1059 * t1112 + t1062 * t1101) * t1074 + (-t1030 * t1059 + t1062 * t1146) * t1088 + t1021 * t1152 + t1059 * t1119) * t1104 - ((t1058 * t1111 + t1061 * t1102) * t1073 + (-t1029 * t1058 + t1061 * t1147) * t1086 + t1020 * t1153 + t1058 * t1120) * t1105 - ((t1057 * t1110 + t1060 * t1103) * t1072 + (-t1028 * t1057 + t1060 * t1148) * t1084 + t1019 * t1154 + t1057 * t1121) * t1106) * m(3); -g(2) * m(4) + (t1059 * t1128 + ((-t1059 * t1125 - t1152) * t1074 + (-t1027 * t1059 + t1062 * t1129) * t1088 + t1062 * t1046) * t1149) * t1033 + (t1058 * t1131 + ((-t1058 * t1126 - t1153) * t1073 + (-t1025 * t1058 + t1061 * t1132) * t1086 + t1061 * t1045) * t1150) * t1032 + (t1057 * t1134 + ((-t1057 * t1127 - t1154) * t1072 + (-t1023 * t1057 + t1060 * t1135) * t1084 + t1060 * t1044) * t1151) * t1031 + (-((-t1059 * t1101 + t1062 * t1112) * t1074 + (-t1030 * t1062 - t1059 * t1146) * t1088 - t1021 * t1155 + t1062 * t1119) * t1104 - ((-t1058 * t1102 + t1061 * t1111) * t1073 + (-t1029 * t1061 - t1058 * t1147) * t1086 - t1020 * t1156 + t1061 * t1120) * t1105 - ((-t1057 * t1103 + t1060 * t1110) * t1072 + (-t1028 * t1060 - t1057 * t1148) * t1084 - t1019 * t1157 + t1060 * t1121) * t1106) * m(3); -t1089 * t1033 * t1012 - t1088 * (t1083 * t1098 - t1122) * t1015 * t1143 - t1087 * t1032 * t1011 - t1086 * (t1081 * t1099 - t1123) * t1014 * t1144 - t1085 * t1031 * t1010 - t1084 * (t1079 * t1100 - t1124) * t1013 * t1145 - g(3) * m(4) + (-(-t1083 * t1074 * t1137 - ((0.2e1 * t1056 + pkin(1)) * t1083 - t1122) * t1116 - qJ(3,1) * (t1046 * t1083 - t1082 * t1122)) * t1104 - (-t1081 * t1073 * t1138 - ((0.2e1 * t1055 + pkin(1)) * t1081 - t1123) * t1117 - qJ(3,2) * (t1045 * t1081 - t1080 * t1123)) * t1105 - (-t1079 * t1072 * t1139 - ((0.2e1 * t1054 + pkin(1)) * t1079 - t1124) * t1118 - qJ(3,3) * (t1044 * t1079 - t1078 * t1124)) * t1106) * m(3);];
taugX  = t1;
