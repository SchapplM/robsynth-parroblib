% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRRRR1V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:33:57
% EndTime: 2020-08-07 03:34:06
% DurationCPUTime: 9.25s
% Computational Cost: add. (40950->546), mult. (73983->886), div. (18090->8), fcn. (80607->36), ass. (0->373)
t935 = cos(qJ(1,3));
t1144 = pkin(2) * mrSges(3,1);
t1002 = t1144 / 0.2e1;
t947 = xDP(1);
t955 = 0.1e1 / pkin(2);
t1055 = t947 * t955;
t946 = xDP(2);
t1056 = t946 * t955;
t916 = qJ(2,3) + qJ(3,3);
t1161 = cos(qJ(1,3) - t916) + cos(qJ(1,3) + t916);
t924 = sin(qJ(3,3));
t904 = 0.1e1 / t924;
t1089 = t1161 * t904;
t921 = legFrame(3,2);
t889 = sin(t921);
t926 = sin(qJ(1,3));
t1086 = t889 * t926;
t892 = cos(t921);
t925 = sin(qJ(2,3));
t1071 = t924 * t925;
t933 = cos(qJ(3,3));
t934 = cos(qJ(2,3));
t974 = -t933 * t934 + t1071;
t1070 = t924 * t934;
t975 = -t925 * t933 - t1070;
t818 = t974 * t1086 + t975 * t892;
t1083 = t892 * t926;
t819 = -t974 * t1083 + t975 * t889;
t945 = xDP(3);
t988 = t945 * t955 / 0.2e1;
t788 = t988 * t1089 + (t1055 * t819 + t1056 * t818) * t904;
t1124 = pkin(3) * t933;
t882 = pkin(2) + t1124;
t839 = pkin(3) * t1070 + t882 * t925;
t1012 = pkin(3) * t1071;
t857 = t882 * t934;
t842 = t857 - t1012;
t824 = -t842 * t1083 + t839 * t889;
t953 = 0.1e1 / pkin(3);
t1054 = t953 * t955;
t992 = t947 * t1054;
t803 = t824 * t904 * t992;
t827 = t842 * t1086 + t839 * t892;
t993 = t946 * t1054;
t806 = t827 * t904 * t993;
t994 = t945 * t1054;
t997 = t842 * t904 * t935;
t968 = t994 * t997;
t773 = t803 / 0.2e1 + t806 / 0.2e1 - t968 / 0.2e1 + t788;
t1139 = pkin(2) * t773;
t1018 = t924 * t1139;
t779 = t803 + t806 - t968;
t776 = t779 + t788;
t907 = t933 ^ 2;
t1106 = t776 * t907;
t1027 = -0.4e1 * t1106;
t1111 = t924 * mrSges(3,1);
t1030 = pkin(2) * t1111;
t1151 = -2 * pkin(1);
t895 = m(3) * pkin(2) + mrSges(2,1);
t1043 = t895 * t1151;
t1046 = mrSges(2,2) * t1151;
t919 = Ifges(3,2) - Ifges(3,1);
t1077 = t919 * t907;
t1105 = t776 * t924;
t1108 = pkin(2) * mrSges(3,3) - Ifges(2,5);
t1114 = t907 * Ifges(3,4);
t1117 = Ifges(3,4) * t924;
t920 = Ifges(2,4) - Ifges(3,4);
t1121 = 0.2e1 * t920;
t1133 = pkin(2) * t788;
t1142 = pkin(1) * t925;
t954 = pkin(2) ^ 2;
t950 = m(3) * t954;
t869 = -Ifges(2,1) + Ifges(2,2) + t950 - t919;
t1146 = -0.2e1 * t869;
t845 = 0.1e1 / (pkin(3) * cos(t916) + t934 * pkin(2) + pkin(1));
t800 = (-t926 * t945 + (-t889 * t946 + t892 * t947) * t935) * t845;
t1149 = -0.2e1 * t800;
t1150 = 2 * pkin(1);
t1158 = 4 * mrSges(3,2);
t1000 = t1133 / 0.2e1;
t897 = t933 * pkin(2);
t1006 = t788 * t897;
t1009 = pkin(3) * t1105;
t1015 = t776 * t1124;
t1068 = t926 * t935;
t1069 = t925 * t934;
t1080 = t904 * t955;
t1107 = t776 * t779;
t1130 = pkin(3) * t776;
t1152 = -0.2e1 * t1130;
t1157 = 0.2e1 * t935 ^ 2;
t761 = t1015 + t1133;
t1127 = pkin(3) * t907;
t854 = -pkin(3) + t897 + 0.2e1 * t1127;
t863 = pkin(1) - 0.2e1 * t1012;
t908 = t934 ^ 2;
t962 = pkin(1) * t1071 - pkin(3) + t1127;
t737 = (-pkin(3) * t1107 + (-0.4e1 * ((t933 * t1000 + (t907 - 0.1e1 / 0.2e1) * t1130) * t1069 + ((t1000 + t1015) * t908 - t761 / 0.2e1) * t924) * t1068 + t1161 * ((t934 * t1009 + t761 * t925) * t926 - t935 * (pkin(1) + t842) * t800) + (t1157 - 0.2e1) * (t854 * t908 + (-pkin(2) * t1071 + t863 * t933) * t934 - t962) * t800) * t800 / 0.2e1 + (((t1006 + (0.2e1 * t907 - 0.1e1) * t1130) * t908 - (0.2e1 * t1015 + t1133) * t924 * t1069 + t1130 - t1130 * t907) * t1157 + (t854 * t1069 + ((pkin(2) + 0.2e1 * t1124) * t908 - t882) * t924) * t1068 * t1149 - 0.2e1 * t1006 + t1152 + t1161 * ((t925 * t1009 - t761 * t934) * t935 + t926 * t839 * t800)) * t788 / 0.2e1) * t1080;
t1042 = pkin(3) * t897;
t952 = pkin(3) ^ 2;
t1047 = t952 - t954;
t1145 = -0.2e1 * t952;
t797 = t800 ^ 2;
t743 = ((t897 + pkin(3)) * t1107 + (-((t907 * t1145 - 0.2e1 * t1042 + t1047) * t908 - t863 * t857 + pkin(3) * t962) * t797 + (0.2e1 * t773 * t1042 + t776 * t952 + t788 * t954) * t788) * t953) * t1080;
t746 = (-0.2e1 * t925 * t1133 + sin(t916) * t1152) * t845 * t800;
t943 = pkin(1) * mrSges(3,2);
t1044 = 2 * t943;
t764 = Ifges(3,5) * t776 + t800 * t1044;
t944 = pkin(1) * mrSges(3,1);
t1045 = -0.2e1 * t944;
t767 = Ifges(3,6) * t776 + t800 * t1045;
t1074 = t919 * t924;
t770 = t776 * t1074;
t782 = t920 * t788;
t785 = t788 ^ 2;
t830 = -(-Ifges(3,5) * t933 + Ifges(3,6) * t924 + t1108) * t925 + (Ifges(3,5) * t924 + Ifges(3,6) * t933 + Ifges(2,6)) * t934;
t833 = -(-Ifges(3,5) * t934 + Ifges(3,6) * t925) * t924 + t933 * (Ifges(3,5) * t925 + Ifges(3,6) * t934);
t1120 = mrSges(3,2) * t924;
t860 = t895 - t1120;
t870 = mrSges(2,2) + t1111;
t1143 = pkin(2) * mrSges(3,2);
t901 = -0.2e1 * t1143;
t903 = 0.2e1 * t944;
t1024 = 0.2e1 * t1077;
t1033 = pkin(2) * t1120;
t1039 = 0.4e1 * t1117;
t902 = 0.2e1 * t1144;
t959 = (t902 + t1039) * t933 + t1024 + t869 - 0.2e1 * t1033;
t969 = -(m(2) + m(3)) * (pkin(1) ^ 2) - Ifges(2,1) - Ifges(3,2) - Ifges(1,3);
t987 = t845 * ((-t959 * t908 - (t860 * t1150 + t903 * t933 + (0.4e1 * t1114 + t901 * t933 + 0.2e1 * (-t919 * t933 - t1144) * t924 + t1121) * t925) * t934 + t1077 - 0.2e1 * (-mrSges(3,2) * t1142 - t1117) * t933 + 0.2e1 * t870 * t1142 + t969) * t746 + t830 * t737 + t833 * t743 + 0.4e1 * t800 * ((-mrSges(3,2) * t1139 - t770) * t933 - mrSges(3,1) * t1018 + t782 + (-t779 + 0.2e1 * t1106) * Ifges(3,4)) * t908 + (t785 * t1108 + (-t764 * t933 + t767 * t924) * t776 + (t788 * t1046 + (-0.8e1 * (Ifges(3,4) * t1105 + t773 * t1002) * t933 + t1018 * t1158 + t788 * t1146 + (0.2e1 * t779 + t1027) * t919) * t925) * t800) * t934 + Ifges(3,4) * t800 * t1027 + (t767 * t776 * t925 + (-mrSges(3,2) * t1133 - t770) * t1149) * t933 + (t764 * t1105 + (Ifges(2,6) * t788 + t800 * t1043) * t788) * t925 + (-Ifges(3,4) * t779 - t788 * t1030 + t782) * t1149);
t1164 = t935 * t987;
t938 = cos(qJ(1,2));
t917 = qJ(2,2) + qJ(3,2);
t1160 = cos(qJ(1,2) - t917) + cos(qJ(1,2) + t917);
t927 = sin(qJ(3,2));
t905 = 0.1e1 / t927;
t1088 = t1160 * t905;
t922 = legFrame(2,2);
t890 = sin(t922);
t929 = sin(qJ(1,2));
t1085 = t890 * t929;
t893 = cos(t922);
t928 = sin(qJ(2,2));
t1067 = t927 * t928;
t936 = cos(qJ(3,2));
t937 = cos(qJ(2,2));
t972 = -t936 * t937 + t1067;
t1066 = t927 * t937;
t973 = -t928 * t936 - t1066;
t820 = t972 * t1085 + t973 * t893;
t1082 = t893 * t929;
t821 = -t972 * t1082 + t973 * t890;
t789 = t988 * t1088 + (t1055 * t821 + t1056 * t820) * t905;
t1123 = pkin(3) * t936;
t883 = pkin(2) + t1123;
t840 = pkin(3) * t1066 + t883 * t928;
t1011 = pkin(3) * t1067;
t858 = t883 * t937;
t843 = t858 - t1011;
t825 = -t843 * t1082 + t840 * t890;
t804 = t825 * t905 * t992;
t828 = t843 * t1085 + t840 * t893;
t807 = t828 * t905 * t993;
t996 = t843 * t905 * t938;
t967 = t994 * t996;
t774 = t804 / 0.2e1 + t807 / 0.2e1 - t967 / 0.2e1 + t789;
t1138 = pkin(2) * t774;
t1017 = t927 * t1138;
t780 = t804 + t807 - t967;
t777 = t780 + t789;
t910 = t936 ^ 2;
t1103 = t777 * t910;
t1026 = -0.4e1 * t1103;
t1110 = t927 * mrSges(3,1);
t1029 = pkin(2) * t1110;
t1076 = t919 * t910;
t1102 = t777 * t927;
t1113 = t910 * Ifges(3,4);
t1116 = Ifges(3,4) * t927;
t1132 = pkin(2) * t789;
t1141 = pkin(1) * t928;
t846 = 0.1e1 / (pkin(3) * cos(t917) + t937 * pkin(2) + pkin(1));
t801 = (-t929 * t945 + (-t890 * t946 + t893 * t947) * t938) * t846;
t1148 = -0.2e1 * t801;
t898 = t936 * pkin(2);
t1005 = t789 * t898;
t1008 = pkin(3) * t1102;
t1014 = t777 * t1123;
t1064 = t929 * t938;
t1065 = t928 * t937;
t1079 = t905 * t955;
t1104 = t777 * t780;
t1129 = pkin(3) * t777;
t1153 = -0.2e1 * t1129;
t1156 = 0.2e1 * t938 ^ 2;
t762 = t1014 + t1132;
t1126 = pkin(3) * t910;
t855 = -pkin(3) + t898 + 0.2e1 * t1126;
t864 = pkin(1) - 0.2e1 * t1011;
t911 = t937 ^ 2;
t961 = pkin(1) * t1067 - pkin(3) + t1126;
t999 = t1132 / 0.2e1;
t738 = (-pkin(3) * t1104 + (-0.4e1 * ((t936 * t999 + (t910 - 0.1e1 / 0.2e1) * t1129) * t1065 + ((t999 + t1014) * t911 - t762 / 0.2e1) * t927) * t1064 + t1160 * ((t937 * t1008 + t762 * t928) * t929 - t938 * (pkin(1) + t843) * t801) + (t1156 - 0.2e1) * (t855 * t911 + (-pkin(2) * t1067 + t864 * t936) * t937 - t961) * t801) * t801 / 0.2e1 + (((t1005 + (0.2e1 * t910 - 0.1e1) * t1129) * t911 - (0.2e1 * t1014 + t1132) * t927 * t1065 + t1129 - t1129 * t910) * t1156 + (t855 * t1065 + ((pkin(2) + 0.2e1 * t1123) * t911 - t883) * t927) * t1064 * t1148 - 0.2e1 * t1005 + t1153 + t1160 * ((t928 * t1008 - t762 * t937) * t938 + t929 * t840 * t801)) * t789 / 0.2e1) * t1079;
t1041 = pkin(3) * t898;
t798 = t801 ^ 2;
t744 = ((t898 + pkin(3)) * t1104 + (-((t910 * t1145 - 0.2e1 * t1041 + t1047) * t911 - t864 * t858 + pkin(3) * t961) * t798 + (0.2e1 * t774 * t1041 + t777 * t952 + t789 * t954) * t789) * t953) * t1079;
t747 = (-0.2e1 * t928 * t1132 + sin(t917) * t1153) * t846 * t801;
t765 = Ifges(3,5) * t777 + t801 * t1044;
t768 = Ifges(3,6) * t777 + t801 * t1045;
t1073 = t919 * t927;
t771 = t777 * t1073;
t783 = t920 * t789;
t786 = t789 ^ 2;
t831 = -(-Ifges(3,5) * t936 + Ifges(3,6) * t927 + t1108) * t928 + (Ifges(3,5) * t927 + Ifges(3,6) * t936 + Ifges(2,6)) * t937;
t834 = -(-Ifges(3,5) * t937 + Ifges(3,6) * t928) * t927 + t936 * (Ifges(3,5) * t928 + Ifges(3,6) * t937);
t1119 = mrSges(3,2) * t927;
t861 = t895 - t1119;
t871 = mrSges(2,2) + t1110;
t1023 = 0.2e1 * t1076;
t1032 = pkin(2) * t1119;
t1038 = 0.4e1 * t1116;
t958 = (t902 + t1038) * t936 + t1023 + t869 - 0.2e1 * t1032;
t986 = t846 * ((-t958 * t911 - (t861 * t1150 + t903 * t936 + (0.4e1 * t1113 + t901 * t936 + 0.2e1 * (-t919 * t936 - t1144) * t927 + t1121) * t928) * t937 + t1076 - 0.2e1 * (-mrSges(3,2) * t1141 - t1116) * t936 + 0.2e1 * t871 * t1141 + t969) * t747 + t831 * t738 + t834 * t744 + 0.4e1 * t801 * ((-mrSges(3,2) * t1138 - t771) * t936 - mrSges(3,1) * t1017 + t783 + (-t780 + 0.2e1 * t1103) * Ifges(3,4)) * t911 + (t786 * t1108 + (-t765 * t936 + t768 * t927) * t777 + (t789 * t1046 + (-0.8e1 * (Ifges(3,4) * t1102 + t774 * t1002) * t936 + t1017 * t1158 + t789 * t1146 + (0.2e1 * t780 + t1026) * t919) * t928) * t801) * t937 + Ifges(3,4) * t801 * t1026 + (t768 * t777 * t928 + (-mrSges(3,2) * t1132 - t771) * t1148) * t936 + (t765 * t1102 + (Ifges(2,6) * t789 + t801 * t1043) * t789) * t928 + (-Ifges(3,4) * t780 - t789 * t1029 + t783) * t1148);
t1163 = t938 * t986;
t941 = cos(qJ(1,1));
t918 = qJ(2,1) + qJ(3,1);
t1159 = cos(qJ(1,1) - t918) + cos(qJ(1,1) + t918);
t930 = sin(qJ(3,1));
t906 = 0.1e1 / t930;
t1087 = t1159 * t906;
t923 = legFrame(1,2);
t891 = sin(t923);
t932 = sin(qJ(1,1));
t1084 = t891 * t932;
t894 = cos(t923);
t931 = sin(qJ(2,1));
t1063 = t930 * t931;
t939 = cos(qJ(3,1));
t940 = cos(qJ(2,1));
t970 = -t939 * t940 + t1063;
t1062 = t930 * t940;
t971 = -t931 * t939 - t1062;
t822 = t970 * t1084 + t971 * t894;
t1081 = t894 * t932;
t823 = -t970 * t1081 + t971 * t891;
t790 = t988 * t1087 + (t1055 * t823 + t1056 * t822) * t906;
t1122 = pkin(3) * t939;
t884 = pkin(2) + t1122;
t841 = pkin(3) * t1062 + t884 * t931;
t1010 = pkin(3) * t1063;
t859 = t884 * t940;
t844 = t859 - t1010;
t826 = -t844 * t1081 + t841 * t891;
t805 = t826 * t906 * t992;
t829 = t844 * t1084 + t841 * t894;
t808 = t829 * t906 * t993;
t995 = t844 * t906 * t941;
t966 = t994 * t995;
t775 = t805 / 0.2e1 + t808 / 0.2e1 - t966 / 0.2e1 + t790;
t1137 = pkin(2) * t775;
t1016 = t930 * t1137;
t781 = t805 + t808 - t966;
t778 = t781 + t790;
t913 = t939 ^ 2;
t1100 = t778 * t913;
t1025 = -0.4e1 * t1100;
t1109 = t930 * mrSges(3,1);
t1028 = pkin(2) * t1109;
t1075 = t919 * t913;
t1099 = t778 * t930;
t1112 = t913 * Ifges(3,4);
t1115 = Ifges(3,4) * t930;
t1131 = pkin(2) * t790;
t1140 = pkin(1) * t931;
t847 = 0.1e1 / (pkin(3) * cos(t918) + t940 * pkin(2) + pkin(1));
t802 = (-t932 * t945 + (-t891 * t946 + t894 * t947) * t941) * t847;
t1147 = -0.2e1 * t802;
t899 = t939 * pkin(2);
t1004 = t790 * t899;
t1007 = pkin(3) * t1099;
t1013 = t778 * t1122;
t1060 = t932 * t941;
t1061 = t931 * t940;
t1078 = t906 * t955;
t1101 = t778 * t781;
t1128 = pkin(3) * t778;
t1154 = -0.2e1 * t1128;
t1155 = 0.2e1 * t941 ^ 2;
t763 = t1013 + t1131;
t1125 = pkin(3) * t913;
t856 = -pkin(3) + t899 + 0.2e1 * t1125;
t865 = pkin(1) - 0.2e1 * t1010;
t914 = t940 ^ 2;
t960 = pkin(1) * t1063 - pkin(3) + t1125;
t998 = t1131 / 0.2e1;
t739 = (-pkin(3) * t1101 + (-0.4e1 * ((t939 * t998 + (t913 - 0.1e1 / 0.2e1) * t1128) * t1061 + ((t998 + t1013) * t914 - t763 / 0.2e1) * t930) * t1060 + t1159 * ((t940 * t1007 + t763 * t931) * t932 - t941 * (pkin(1) + t844) * t802) + (t1155 - 0.2e1) * (t856 * t914 + (-pkin(2) * t1063 + t865 * t939) * t940 - t960) * t802) * t802 / 0.2e1 + (((t1004 + (0.2e1 * t913 - 0.1e1) * t1128) * t914 - (0.2e1 * t1013 + t1131) * t930 * t1061 + t1128 - t1128 * t913) * t1155 + (t856 * t1061 + ((pkin(2) + 0.2e1 * t1122) * t914 - t884) * t930) * t1060 * t1147 - 0.2e1 * t1004 + t1154 + t1159 * ((t931 * t1007 - t763 * t940) * t941 + t932 * t841 * t802)) * t790 / 0.2e1) * t1078;
t1040 = pkin(3) * t899;
t799 = t802 ^ 2;
t745 = ((t899 + pkin(3)) * t1101 + (-((t913 * t1145 - 0.2e1 * t1040 + t1047) * t914 - t865 * t859 + pkin(3) * t960) * t799 + (0.2e1 * t775 * t1040 + t778 * t952 + t790 * t954) * t790) * t953) * t1078;
t748 = (-0.2e1 * t931 * t1131 + sin(t918) * t1154) * t847 * t802;
t766 = Ifges(3,5) * t778 + t802 * t1044;
t769 = Ifges(3,6) * t778 + t802 * t1045;
t1072 = t919 * t930;
t772 = t778 * t1072;
t784 = t920 * t790;
t787 = t790 ^ 2;
t832 = -(-Ifges(3,5) * t939 + Ifges(3,6) * t930 + t1108) * t931 + (Ifges(3,5) * t930 + Ifges(3,6) * t939 + Ifges(2,6)) * t940;
t835 = -(-Ifges(3,5) * t940 + Ifges(3,6) * t931) * t930 + t939 * (Ifges(3,5) * t931 + Ifges(3,6) * t940);
t1118 = mrSges(3,2) * t930;
t862 = t895 - t1118;
t872 = mrSges(2,2) + t1109;
t1022 = 0.2e1 * t1075;
t1031 = pkin(2) * t1118;
t1037 = 0.4e1 * t1115;
t957 = (t902 + t1037) * t939 + t1022 + t869 - 0.2e1 * t1031;
t985 = t847 * ((-t957 * t914 - (t862 * t1150 + t903 * t939 + (0.4e1 * t1112 + t901 * t939 + 0.2e1 * (-t919 * t939 - t1144) * t930 + t1121) * t931) * t940 + t1075 - 0.2e1 * (-mrSges(3,2) * t1140 - t1115) * t939 + 0.2e1 * t872 * t1140 + t969) * t748 + t832 * t739 + t835 * t745 + 0.4e1 * t802 * ((-mrSges(3,2) * t1137 - t772) * t939 - mrSges(3,1) * t1016 + t784 + (-t781 + 0.2e1 * t1100) * Ifges(3,4)) * t914 + (t787 * t1108 + (-t766 * t939 + t769 * t930) * t778 + (t790 * t1046 + (-0.8e1 * (Ifges(3,4) * t1099 + t775 * t1002) * t939 + t1016 * t1158 + t790 * t1146 + (0.2e1 * t781 + t1025) * t919) * t931) * t802) * t940 + Ifges(3,4) * t802 * t1025 + (t769 * t778 * t931 + (-mrSges(3,2) * t1131 - t772) * t1147) * t939 + (t766 * t1099 + (Ifges(2,6) * t790 + t802 * t1043) * t790) * t931 + (-Ifges(3,4) * t781 - t790 * t1028 + t784) * t1147);
t1162 = t941 * t985;
t1136 = pkin(2) * t785;
t1135 = pkin(2) * t786;
t1134 = pkin(2) * t787;
t1098 = t797 * t908;
t1097 = t797 * t934;
t1096 = t798 * t911;
t1095 = t798 * t937;
t1094 = t799 * t914;
t1093 = t799 * t940;
t1059 = t943 * t933;
t1058 = t943 * t936;
t1057 = t943 * t939;
t1003 = -t1144 / 0.4e1;
t1021 = t797 * t1142;
t1036 = 0.2e1 * t1114;
t794 = t797 * t1036;
t965 = (-mrSges(3,1) * t933 + t1120) * pkin(2);
t848 = -Ifges(3,3) + t965;
t896 = Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1;
t900 = -t1143 / 0.4e1;
t949 = -Ifges(3,4) / 0.2e1;
t978 = mrSges(3,1) * t1021 - t797 * t1074;
t1050 = -Ifges(3,3) * t743 + t737 * t848 + t746 * t833 - 0.4e1 * (t1114 + (t896 * t924 + t900) * t933 + t924 * t1003 + t949) * t1098 + (t944 * t924 + t1059 + (t1024 + (t1039 + t1144) * t933 - t1033 - t919) * t925) * t1097 + t794 + (mrSges(3,2) * t1136 + t978) * t933 + (mrSges(3,1) * t1136 - mrSges(3,2) * t1021) * t924 - Ifges(3,4) * t797;
t1020 = t798 * t1141;
t1035 = 0.2e1 * t1113;
t795 = t798 * t1035;
t964 = (-mrSges(3,1) * t936 + t1119) * pkin(2);
t849 = -Ifges(3,3) + t964;
t977 = mrSges(3,1) * t1020 - t798 * t1073;
t1049 = -Ifges(3,3) * t744 + t738 * t849 + t747 * t834 - 0.4e1 * (t1113 + (t896 * t927 + t900) * t936 + t927 * t1003 + t949) * t1096 + (t944 * t927 + t1058 + (t1023 + (t1038 + t1144) * t936 - t1032 - t919) * t928) * t1095 + t795 + (mrSges(3,2) * t1135 + t977) * t936 + (mrSges(3,1) * t1135 - mrSges(3,2) * t1020) * t927 - Ifges(3,4) * t798;
t1019 = t799 * t1140;
t1034 = 0.2e1 * t1112;
t796 = t799 * t1034;
t963 = (-mrSges(3,1) * t939 + t1118) * pkin(2);
t850 = -Ifges(3,3) + t963;
t976 = mrSges(3,1) * t1019 - t799 * t1072;
t1048 = -Ifges(3,3) * t745 + t739 * t850 + t748 * t835 - 0.4e1 * (t1112 + (t896 * t930 + t900) * t939 + t930 * t1003 + t949) * t1094 + (t944 * t930 + t1057 + (t1022 + (t1037 + t1144) * t939 - t1031 - t919) * t931) * t1093 + t796 + (mrSges(3,2) * t1134 + t976) * t939 + (mrSges(3,1) * t1134 - mrSges(3,2) * t1019) * t930 - Ifges(3,4) * t799;
t1001 = -Ifges(2,3) - Ifges(3,3) - t950;
t728 = t830 * t746 + (t1001 + 0.2e1 * t965) * t737 + t848 * t743;
t758 = t797 + (0.2e1 * t788 + t779) * t779;
t740 = -0.2e1 * (t1036 + (-t1074 - t1143) * t933 - t1030 + t920) * t1098 + (pkin(1) * t870 + t959 * t925 + t1059) * t1097 + t794 + (-t758 * t1143 + t978) * t933 + t860 * t1021 - t758 * t1030 + t797 * t920;
t984 = (t728 + t740) * t904;
t729 = t831 * t747 + (t1001 + 0.2e1 * t964) * t738 + t849 * t744;
t759 = t798 + (0.2e1 * t789 + t780) * t780;
t741 = -0.2e1 * (t1035 + (-t1073 - t1143) * t936 - t1029 + t920) * t1096 + (pkin(1) * t871 + t958 * t928 + t1058) * t1095 + t795 + (-t759 * t1143 + t977) * t936 + t861 * t1020 - t759 * t1029 + t798 * t920;
t983 = (t729 + t741) * t905;
t730 = t832 * t748 + (t1001 + 0.2e1 * t963) * t739 + t850 * t745;
t760 = t799 + (0.2e1 * t790 + t781) * t781;
t742 = -0.2e1 * (t1034 + (-t1072 - t1143) * t939 - t1028 + t920) * t1094 + (pkin(1) * t872 + t957 * t931 + t1057) * t1093 + t796 + (-t760 * t1143 + t976) * t939 + t862 * t1019 - t760 * t1028 + t799 * t920;
t982 = (t730 + t742) * t906;
t981 = t1050 * t904;
t980 = t1049 * t905;
t979 = t1048 * t906;
t1 = [t894 * t1162 + t893 * t1163 + t892 * t1164 + (t823 * t982 + t821 * t983 + t819 * t984 + (t824 * t981 + t825 * t980 + t826 * t979) * t953) * t955; -t891 * t1162 - t890 * t1163 - t889 * t1164 + (t822 * t982 + t820 * t983 + t818 * t984 + (t827 * t981 + t828 * t980 + t829 * t979) * t953) * t955; -t932 * t985 - t929 * t986 - t926 * t987 + ((t730 / 0.2e1 + t742 / 0.2e1) * t1087 + (t729 / 0.2e1 + t741 / 0.2e1) * t1088 + (t728 / 0.2e1 + t740 / 0.2e1) * t1089 + (-t1048 * t995 - t1049 * t996 - t1050 * t997) * t953) * t955;];
taucX  = t1;
