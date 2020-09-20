% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRRRR2G2P2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2P2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:08:38
% EndTime: 2020-03-09 21:08:41
% DurationCPUTime: 3.46s
% Computational Cost: add. (12321->285), mult. (24504->538), div. (9279->14), fcn. (25965->51), ass. (0->274)
t893 = sin(qJ(1,1));
t1033 = pkin(1) * t893;
t884 = legFrame(1,2);
t854 = sin(t884);
t857 = cos(t884);
t900 = cos(qJ(3,1));
t892 = sin(qJ(2,1));
t901 = cos(qJ(2,1));
t902 = cos(qJ(1,1));
t820 = t892 * t902 + t893 * t901;
t873 = t900 ^ 2;
t970 = pkin(2) * t820 * t873;
t1030 = pkin(1) * t901;
t891 = sin(qJ(3,1));
t973 = t891 * t1030;
t994 = t857 * t891;
t777 = t854 * t970 + (-pkin(2) * t994 + t854 * t1033) * t900 - t857 * t973;
t908 = xDP(2);
t875 = 0.1e1 / t900 ^ 2;
t866 = 0.1e1 / t892;
t917 = 0.1e1 / pkin(2);
t919 = 0.1e1 / pkin(1);
t983 = t917 * t919;
t955 = t866 * t983;
t941 = t875 * t955;
t771 = t777 * t908 * t941;
t997 = t854 * t891;
t778 = -t857 * t970 + (-pkin(2) * t997 - t857 * t1033) * t900 - t854 * t973;
t909 = xDP(1);
t772 = t778 * t909 * t941;
t1049 = -0.2e1 * pkin(1);
t880 = qJ(2,1) + qJ(3,1);
t841 = sin(t880);
t881 = qJ(2,1) - qJ(3,1);
t842 = sin(t881);
t1009 = (t902 * t1049 + (-cos(qJ(1,1) + t881) - cos(qJ(1,1) + t880)) * pkin(2)) / (t841 + t842);
t907 = xDP(3);
t954 = t907 * t983;
t787 = t954 * t1009;
t757 = t772 + t771 + t787;
t1006 = t820 * t900;
t795 = -t854 * t1006 + t994;
t796 = t857 * t1006 + t997;
t851 = cos(qJ(1,1) + qJ(2,1));
t874 = 0.1e1 / t900;
t991 = t866 * t919;
t956 = t874 * t991;
t984 = t907 * t919;
t763 = t851 * t866 * t984 + (t795 * t908 + t796 * t909) * t956;
t977 = -t757 - t763;
t1056 = t977 ^ 2;
t890 = sin(qJ(1,2));
t1034 = pkin(1) * t890;
t883 = legFrame(2,2);
t853 = sin(t883);
t856 = cos(t883);
t897 = cos(qJ(3,2));
t889 = sin(qJ(2,2));
t898 = cos(qJ(2,2));
t899 = cos(qJ(1,2));
t819 = t889 * t899 + t890 * t898;
t870 = t897 ^ 2;
t971 = pkin(2) * t819 * t870;
t1031 = pkin(1) * t898;
t888 = sin(qJ(3,2));
t974 = t888 * t1031;
t995 = t856 * t888;
t775 = t853 * t971 + (-pkin(2) * t995 + t853 * t1034) * t897 - t856 * t974;
t872 = 0.1e1 / t897 ^ 2;
t865 = 0.1e1 / t889;
t957 = t865 * t983;
t942 = t872 * t957;
t769 = t775 * t908 * t942;
t998 = t853 * t888;
t776 = -t856 * t971 + (-pkin(2) * t998 - t856 * t1034) * t897 - t853 * t974;
t770 = t776 * t909 * t942;
t878 = qJ(2,2) + qJ(3,2);
t839 = sin(t878);
t879 = qJ(2,2) - qJ(3,2);
t840 = sin(t879);
t1010 = (t899 * t1049 + (-cos(qJ(1,2) + t879) - cos(qJ(1,2) + t878)) * pkin(2)) / (t839 + t840);
t786 = t954 * t1010;
t756 = t770 + t769 + t786;
t1007 = t819 * t897;
t793 = -t853 * t1007 + t995;
t794 = t856 * t1007 + t998;
t848 = cos(qJ(1,2) + qJ(2,2));
t871 = 0.1e1 / t897;
t992 = t865 * t919;
t958 = t871 * t992;
t762 = t848 * t865 * t984 + (t793 * t908 + t794 * t909) * t958;
t978 = -t756 - t762;
t1055 = t978 ^ 2;
t887 = sin(qJ(1,3));
t1035 = pkin(1) * t887;
t882 = legFrame(3,2);
t852 = sin(t882);
t855 = cos(t882);
t894 = cos(qJ(3,3));
t886 = sin(qJ(2,3));
t895 = cos(qJ(2,3));
t896 = cos(qJ(1,3));
t818 = t886 * t896 + t887 * t895;
t867 = t894 ^ 2;
t972 = pkin(2) * t818 * t867;
t1032 = pkin(1) * t895;
t885 = sin(qJ(3,3));
t975 = t885 * t1032;
t996 = t855 * t885;
t773 = t852 * t972 + (-pkin(2) * t996 + t852 * t1035) * t894 - t855 * t975;
t869 = 0.1e1 / t894 ^ 2;
t864 = 0.1e1 / t886;
t959 = t864 * t983;
t943 = t869 * t959;
t767 = t773 * t908 * t943;
t999 = t852 * t885;
t774 = -t855 * t972 + (-pkin(2) * t999 - t855 * t1035) * t894 - t852 * t975;
t768 = t774 * t909 * t943;
t876 = qJ(2,3) + qJ(3,3);
t837 = sin(t876);
t877 = qJ(2,3) - qJ(3,3);
t838 = sin(t877);
t1011 = (t896 * t1049 + (-cos(qJ(1,3) + t877) - cos(qJ(1,3) + t876)) * pkin(2)) / (t837 + t838);
t785 = t954 * t1011;
t755 = t768 + t767 + t785;
t1008 = t818 * t894;
t791 = -t852 * t1008 + t996;
t792 = t855 * t1008 + t999;
t845 = cos(qJ(1,3) + qJ(2,3));
t868 = 0.1e1 / t894;
t993 = t864 * t919;
t960 = t868 * t993;
t761 = t845 * t864 * t984 + (t791 * t908 + t792 * t909) * t960;
t979 = -t755 - t761;
t1054 = t979 ^ 2;
t1053 = 0.2e1 * pkin(1);
t914 = rSges(3,2) ^ 2;
t915 = rSges(3,1) ^ 2;
t830 = (-t914 + t915) * m(3) - Icges(3,1) + Icges(3,2);
t1045 = -t830 / 0.2e1;
t836 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t911 = 0.2e1 * qJ(3,3);
t858 = sin(t911);
t861 = cos(t911);
t1052 = t861 * t1045 + t836 * t858;
t912 = 0.2e1 * qJ(3,2);
t859 = sin(t912);
t862 = cos(t912);
t1051 = t862 * t1045 + t836 * t859;
t913 = 0.2e1 * qJ(3,1);
t860 = sin(t913);
t863 = cos(t913);
t1050 = t863 * t1045 + t836 * t860;
t758 = t761 ^ 2;
t759 = t762 ^ 2;
t760 = t763 ^ 2;
t1048 = -0.2e1 * t830;
t1047 = -m(3) / 0.2e1;
t906 = m(2) * rSges(2,1);
t1046 = m(3) * rSges(3,3);
t1044 = t837 / 0.2e1;
t1043 = t839 / 0.2e1;
t1042 = t841 / 0.2e1;
t1041 = cos(t876) / 0.2e1;
t1040 = cos(t878) / 0.2e1;
t1039 = cos(t880) / 0.2e1;
t1038 = pkin(1) * t758;
t1037 = pkin(1) * t759;
t1036 = pkin(1) * t760;
t746 = t768 / 0.2e1 + t767 / 0.2e1 + t785 / 0.2e1 + t761;
t1020 = t755 * t746;
t800 = (t852 * t909 + t855 * t908) * t917 * t868;
t797 = t800 ^ 2;
t1029 = (t797 / 0.2e1 + t1020) * t886;
t747 = t770 / 0.2e1 + t769 / 0.2e1 + t786 / 0.2e1 + t762;
t1019 = t756 * t747;
t801 = (t853 * t909 + t856 * t908) * t917 * t871;
t798 = t801 ^ 2;
t1028 = (t798 / 0.2e1 + t1019) * t889;
t748 = t772 / 0.2e1 + t771 / 0.2e1 + t787 / 0.2e1 + t763;
t1018 = t757 * t748;
t802 = (t854 * t909 + t857 * t908) * t917 * t874;
t799 = t802 ^ 2;
t1027 = (t799 / 0.2e1 + t1018) * t892;
t1026 = t979 * t800;
t1025 = t979 * t867;
t1024 = t978 * t801;
t1023 = t978 * t870;
t1022 = t977 * t802;
t1021 = t977 * t873;
t1017 = t797 * t868;
t1016 = t798 * t871;
t1015 = t799 * t874;
t1014 = t800 * t895;
t1013 = t801 * t898;
t1012 = t802 * t901;
t1005 = t830 * t858;
t1004 = t830 * t859;
t1003 = t830 * t860;
t1002 = t836 * t861;
t1001 = t836 * t862;
t1000 = t836 * t863;
t990 = t885 * t886;
t989 = t888 * t889;
t988 = t891 * t892;
t987 = t894 * t895;
t986 = t897 * t898;
t985 = t900 * t901;
t916 = pkin(2) ^ 2;
t963 = t800 * t990;
t731 = (-t916 * t1025 + (pkin(1) * t761 + (0.2e1 * t746 * t987 - t963) * pkin(2)) * pkin(1)) * t868 * t761 * t959 + (-pkin(2) * t1025 + (-t979 * t987 - t963) * pkin(1)) * t755 * t960 + ((pkin(1) * t979 * t990 + pkin(2) * t800) * t894 + pkin(1) * t1014) * t869 * t800 * t993;
t737 = (-t758 * t1032 + (-t894 * t1054 - t1017) * pkin(2)) * t993;
t940 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t976 = t914 + t915;
t953 = 0.2e1 * rSges(3,3) ^ 2 + t976;
t933 = t953 * t1047 + t940;
t923 = t933 + t1052;
t833 = m(2) * rSges(2,2) - t1046;
t929 = ((-t906 + (-rSges(3,1) * t894 + rSges(3,2) * t885) * m(3)) * t895 + t833 * t886) * pkin(1);
t764 = t929 + t923;
t834 = -rSges(3,2) * t1046 + Icges(3,6);
t835 = rSges(3,1) * t1046 - Icges(3,5);
t806 = -t834 * t894 + t835 * t885;
t844 = cos(t877);
t926 = (-t834 * t885 - t835 * t894) * t797;
t932 = -t833 * t895 - t886 * t906;
t966 = t885 * t1017;
t982 = t923 * t731 + t764 * t737 - t806 * t966 + t926 - (-0.2e1 * t1002 - t1005) * t1026 + (((-t844 / 0.2e1 + t1041) * rSges(3,2) + (t838 / 0.2e1 + t1044) * rSges(3,1)) * m(3) - t932) * t1038;
t962 = t801 * t989;
t732 = (-t916 * t1023 + (pkin(1) * t762 + (0.2e1 * t747 * t986 - t962) * pkin(2)) * pkin(1)) * t871 * t762 * t957 + (-pkin(2) * t1023 + (-t978 * t986 - t962) * pkin(1)) * t756 * t958 + ((pkin(1) * t978 * t989 + pkin(2) * t801) * t897 + pkin(1) * t1013) * t872 * t801 * t992;
t738 = (-t759 * t1031 + (-t897 * t1055 - t1016) * pkin(2)) * t992;
t922 = t933 + t1051;
t928 = ((-t906 + (-rSges(3,1) * t897 + rSges(3,2) * t888) * m(3)) * t898 + t833 * t889) * pkin(1);
t765 = t928 + t922;
t807 = -t834 * t897 + t835 * t888;
t847 = cos(t879);
t925 = (-t834 * t888 - t835 * t897) * t798;
t931 = -t833 * t898 - t889 * t906;
t965 = t888 * t1016;
t981 = t922 * t732 + t765 * t738 - t807 * t965 + t925 - (-0.2e1 * t1001 - t1004) * t1024 + (((-t847 / 0.2e1 + t1040) * rSges(3,2) + (t840 / 0.2e1 + t1043) * rSges(3,1)) * m(3) - t931) * t1037;
t961 = t802 * t988;
t733 = (-t916 * t1021 + (pkin(1) * t763 + (0.2e1 * t748 * t985 - t961) * pkin(2)) * pkin(1)) * t874 * t763 * t955 + (-pkin(2) * t1021 + (-t977 * t985 - t961) * pkin(1)) * t757 * t956 + ((pkin(1) * t977 * t988 + pkin(2) * t802) * t900 + pkin(1) * t1012) * t875 * t802 * t991;
t739 = (-t760 * t1030 + (-t900 * t1056 - t1015) * pkin(2)) * t991;
t921 = t933 + t1050;
t927 = ((-t906 + (-rSges(3,1) * t900 + rSges(3,2) * t891) * m(3)) * t901 + t833 * t892) * pkin(1);
t766 = t927 + t921;
t808 = -t834 * t900 + t835 * t891;
t850 = cos(t881);
t924 = (-t834 * t891 - t835 * t900) * t799;
t930 = -t833 * t901 - t892 * t906;
t964 = t891 * t1015;
t980 = t921 * t733 + t766 * t739 - t808 * t964 + t924 - (-0.2e1 * t1000 - t1003) * t1022 + (((-t850 / 0.2e1 + t1039) * rSges(3,2) + (t842 / 0.2e1 + t1042) * rSges(3,1)) * m(3) - t930) * t1036;
t969 = t979 * t1014;
t968 = t978 * t1013;
t967 = t977 * t1012;
t788 = -t885 * Icges(3,5) - Icges(3,6) * t894 + (rSges(3,1) * t885 + rSges(3,2) * t894) * m(3) * (pkin(1) * t886 + rSges(3,3));
t918 = pkin(1) ^ 2;
t920 = -m(2) * t918 - Icges(1,3) + (0.2e1 * t918 + t953) * t1047 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t940;
t949 = (t926 - (t885 * t894 * t1048 + (-0.4e1 * t867 + 0.2e1) * t836) * t1026 + (t932 * t1020 + ((-rSges(3,1) * t1029 + rSges(3,2) * t969) * t894 + (rSges(3,1) * t969 + rSges(3,2) * t1029) * t885) * m(3)) * t1053 + (t920 + 0.2e1 * t929 + t1052) * t737 + t764 * t731 - t788 * t966) * t864;
t789 = -t888 * Icges(3,5) - Icges(3,6) * t897 + (rSges(3,1) * t888 + rSges(3,2) * t897) * m(3) * (pkin(1) * t889 + rSges(3,3));
t948 = (t925 - (t888 * t897 * t1048 + (-0.4e1 * t870 + 0.2e1) * t836) * t1024 + (t931 * t1019 + ((-rSges(3,1) * t1028 + rSges(3,2) * t968) * t897 + (rSges(3,1) * t968 + rSges(3,2) * t1028) * t888) * m(3)) * t1053 + (t920 + 0.2e1 * t928 + t1051) * t738 + t765 * t732 - t789 * t965) * t865;
t790 = -t891 * Icges(3,5) - Icges(3,6) * t900 + (rSges(3,1) * t891 + rSges(3,2) * t900) * m(3) * (pkin(1) * t892 + rSges(3,3));
t947 = (t924 - (t891 * t900 * t1048 + (-0.4e1 * t873 + 0.2e1) * t836) * t1022 + (t930 * t1018 + ((-rSges(3,1) * t1027 + rSges(3,2) * t967) * t900 + (rSges(3,1) * t967 + rSges(3,2) * t1027) * t891) * m(3)) * t1053 + (t920 + 0.2e1 * t927 + t1050) * t739 + t766 * t733 - t790 * t964) * t866;
t832 = -t976 * m(3) - Icges(3,3);
t946 = (t731 * t806 + t737 * t788 - t832 * t966 + (t1005 / 0.2e1 + t1002) * t1054 + ((t844 / 0.2e1 + t1041) * rSges(3,2) + (-t838 / 0.2e1 + t1044) * rSges(3,1)) * m(3) * t1038) * t868;
t945 = (t732 * t807 + t738 * t789 - t832 * t965 + (t1004 / 0.2e1 + t1001) * t1055 + ((t847 / 0.2e1 + t1040) * rSges(3,2) + (-t840 / 0.2e1 + t1043) * rSges(3,1)) * m(3) * t1037) * t871;
t944 = (t733 * t808 + t739 * t790 - t832 * t964 + (t1003 / 0.2e1 + t1000) * t1056 + ((t850 / 0.2e1 + t1039) * rSges(3,2) + (-t842 / 0.2e1 + t1042) * rSges(3,1)) * m(3) * t1036) * t874;
t939 = t868 * t949;
t938 = t871 * t948;
t937 = t874 * t947;
t936 = t982 * t869 * t864;
t935 = t981 * t872 * t865;
t934 = t980 * t875 * t866;
t1 = [(t792 * t939 + t794 * t938 + t796 * t937) * t919 + (t854 * t944 + t853 * t945 + t852 * t946 + (t774 * t936 + t776 * t935 + t778 * t934) * t919) * t917; (t791 * t939 + t793 * t938 + t795 * t937) * t919 + (t857 * t944 + t856 * t945 + t855 * t946 + (t773 * t936 + t775 * t935 + t777 * t934) * t919) * t917; (t851 * t947 + t848 * t948 + t845 * t949 + (t980 * t1009 + t981 * t1010 + t982 * t1011) * t917) * t919;];
taucX  = t1;
