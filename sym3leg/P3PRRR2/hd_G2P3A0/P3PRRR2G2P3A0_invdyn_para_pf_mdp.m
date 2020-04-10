% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRR2G2P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR2G2P3A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR2G2P3A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:42
% EndTime: 2020-03-09 21:21:44
% DurationCPUTime: 2.06s
% Computational Cost: add. (10719->267), mult. (11627->463), div. (2781->10), fcn. (9129->36), ass. (0->241)
t926 = cos(qJ(3,3));
t1032 = pkin(1) * t926;
t936 = 0.1e1 / pkin(2);
t1045 = (pkin(2) + t1032) * t936;
t928 = cos(qJ(3,2));
t1031 = pkin(1) * t928;
t1044 = (pkin(2) + t1031) * t936;
t930 = cos(qJ(3,1));
t1030 = pkin(1) * t930;
t1043 = (pkin(2) + t1030) * t936;
t932 = xDP(3);
t1042 = -2 * t932;
t1041 = -g(1) / 0.2e1;
t1040 = g(2) / 0.2e1;
t911 = qJ(2,3) + qJ(3,3);
t914 = legFrame(3,2);
t884 = t914 + t911;
t869 = sin(t884);
t1039 = t869 / 0.2e1;
t912 = qJ(2,2) + qJ(3,2);
t915 = legFrame(2,2);
t886 = t915 + t912;
t871 = sin(t886);
t1038 = t871 / 0.2e1;
t913 = qJ(2,1) + qJ(3,1);
t916 = legFrame(1,2);
t888 = t916 + t913;
t873 = sin(t888);
t1037 = t873 / 0.2e1;
t885 = -t914 + t911;
t876 = cos(t885);
t1036 = t876 / 0.2e1;
t887 = -t915 + t912;
t878 = cos(t887);
t1035 = t878 / 0.2e1;
t889 = -t916 + t913;
t880 = cos(t889);
t1034 = t880 / 0.2e1;
t938 = 0.1e1 / pkin(1);
t1033 = t938 / 0.4e1;
t1029 = pkin(2) * t926;
t1028 = pkin(2) * t928;
t1027 = pkin(2) * t930;
t1026 = pkin(2) * t932;
t918 = xDDP(2);
t1025 = t918 - g(2);
t919 = xDDP(1);
t1024 = t919 - g(1);
t870 = sin(t885);
t875 = cos(t884);
t896 = sin(t911);
t933 = xDP(2);
t934 = xDP(1);
t809 = t896 * t1042 + (t869 - t870) * t934 + (t875 + t876) * t933;
t920 = sin(qJ(3,3));
t905 = 0.1e1 / t920;
t1023 = t809 * t905;
t872 = sin(t887);
t877 = cos(t886);
t897 = sin(t912);
t810 = t897 * t1042 + (t871 - t872) * t934 + (t877 + t878) * t933;
t922 = sin(qJ(3,2));
t907 = 0.1e1 / t922;
t1022 = t810 * t907;
t874 = sin(t889);
t879 = cos(t888);
t898 = sin(t913);
t811 = t898 * t1042 + (t873 - t874) * t934 + (t879 + t880) * t933;
t924 = sin(qJ(3,1));
t909 = 0.1e1 / t924;
t1021 = t811 * t909;
t921 = sin(qJ(2,3));
t1020 = (pkin(1) * t921 + pkin(2) * t896) * t905;
t923 = sin(qJ(2,2));
t1019 = (pkin(1) * t923 + pkin(2) * t897) * t907;
t925 = sin(qJ(2,1));
t1018 = (pkin(1) * t925 + pkin(2) * t898) * t909;
t1017 = t896 * t905;
t917 = xDDP(3);
t1016 = t896 * t917;
t1015 = t897 * t907;
t1014 = t897 * t917;
t1013 = t898 * t909;
t1012 = t898 * t917;
t899 = sin(t914);
t1011 = t899 * t905;
t1010 = t899 * t919;
t900 = sin(t915);
t1009 = t900 * t907;
t1008 = t900 * t919;
t901 = sin(t916);
t1007 = t901 * t909;
t1006 = t901 * t919;
t902 = cos(t914);
t1005 = t902 * t905;
t1004 = t902 * t918;
t903 = cos(t915);
t1003 = t903 * t907;
t1002 = t903 * t918;
t904 = cos(t916);
t1001 = t904 * t909;
t1000 = t904 * t918;
t999 = t905 * t938;
t906 = 0.1e1 / t920 ^ 2;
t939 = 0.1e1 / pkin(1) ^ 2;
t998 = t906 * t939;
t997 = t907 * t938;
t908 = 0.1e1 / t922 ^ 2;
t996 = t908 * t939;
t995 = t909 * t938;
t910 = 0.1e1 / t924 ^ 2;
t994 = t910 * t939;
t993 = t917 * t938;
t992 = t920 * t921;
t991 = t922 * t923;
t990 = t924 * t925;
t989 = t936 * t938;
t960 = t1023 / 0.2e1;
t803 = t938 * t960;
t833 = t899 * t934 + t902 * t933;
t893 = pkin(1) + t1029;
t927 = cos(qJ(2,3));
t800 = (t920 * t1026 - t833 * t893) * t927 + t921 * (pkin(2) * t833 * t920 + t893 * t932);
t984 = t800 * t989;
t957 = t905 * t984;
t794 = t803 + t957 / 0.2e1;
t797 = t803 + t957;
t935 = pkin(2) ^ 2;
t951 = -t809 * t998 / 0.2e1;
t964 = t917 * t989;
t988 = (t797 * t935 + (0.2e1 * t794 * t1029 + t960) * pkin(1)) * t936 * t951 + t964 * t1020;
t959 = t1022 / 0.2e1;
t804 = t938 * t959;
t834 = t900 * t934 + t903 * t933;
t894 = pkin(1) + t1028;
t929 = cos(qJ(2,2));
t801 = (t922 * t1026 - t834 * t894) * t929 + t923 * (pkin(2) * t834 * t922 + t894 * t932);
t982 = t801 * t989;
t956 = t907 * t982;
t795 = t804 + t956 / 0.2e1;
t798 = t804 + t956;
t950 = -t810 * t996 / 0.2e1;
t987 = (t798 * t935 + (0.2e1 * t795 * t1028 + t959) * pkin(1)) * t936 * t950 + t964 * t1019;
t958 = t1021 / 0.2e1;
t805 = t938 * t958;
t835 = t901 * t934 + t904 * t933;
t895 = pkin(1) + t1027;
t931 = cos(qJ(2,1));
t802 = (t924 * t1026 - t835 * t895) * t931 + t925 * (pkin(2) * t835 * t924 + t895 * t932);
t980 = t802 * t989;
t955 = t909 * t980;
t796 = t805 + t955 / 0.2e1;
t799 = t805 + t955;
t949 = -t811 * t994 / 0.2e1;
t986 = (t799 * t935 + (0.2e1 * t796 * t1027 + t958) * pkin(1)) * t936 * t949 + t964 * t1018;
t985 = t800 * t998;
t983 = t801 * t996;
t981 = t802 * t994;
t830 = -pkin(2) * t992 + t893 * t927;
t979 = t830 * t1011;
t978 = t830 * t1005;
t831 = -pkin(2) * t991 + t894 * t929;
t977 = t831 * t1009;
t976 = t831 * t1003;
t832 = -pkin(2) * t990 + t895 * t931;
t975 = t832 * t1007;
t974 = t832 * t1001;
t836 = t926 * t927 - t992;
t973 = t836 * t1011;
t972 = t836 * t1005;
t971 = t836 * t999;
t837 = t928 * t929 - t991;
t970 = t837 * t1009;
t969 = t837 * t1003;
t968 = t837 * t997;
t838 = t930 * t931 - t990;
t967 = t838 * t1007;
t966 = t838 * t1001;
t965 = t838 * t995;
t963 = t809 ^ 2 * t1033;
t962 = t810 ^ 2 * t1033;
t961 = t811 ^ 2 * t1033;
t788 = t797 * t985;
t789 = t798 * t983;
t790 = t799 * t981;
t791 = -t926 * t1023 / 0.2e1 - pkin(2) * t797;
t815 = t971 * t1010;
t818 = t971 * t1004;
t954 = t791 * t951 + t788 + t815 + t818;
t792 = -t928 * t1022 / 0.2e1 - pkin(2) * t798;
t816 = t968 * t1008;
t819 = t968 * t1002;
t953 = t792 * t950 + t789 + t816 + t819;
t793 = -t930 * t1021 / 0.2e1 - t799 * pkin(2);
t817 = t965 * t1006;
t820 = t965 * t1000;
t952 = t793 * t949 + t790 + t817 + t820;
t948 = g(1) * t1036 + g(2) * t1039 + t870 * t1040 + t875 * t1041 + g(3) * cos(t911);
t947 = g(1) * t1035 + g(2) * t1038 + t872 * t1040 + t877 * t1041 + g(3) * cos(t912);
t946 = g(1) * t1034 + g(2) * t1037 + t874 * t1040 + t879 * t1041 + g(3) * cos(t913);
t945 = g(1) * t1039 + g(2) * t1036 - g(3) * t896 + t875 * t1040 + t870 * t1041;
t944 = g(1) * t1038 + g(2) * t1035 - g(3) * t897 + t877 * t1040 + t872 * t1041;
t943 = g(1) * t1037 + g(2) * t1034 - g(3) * t898 + t879 * t1040 + t874 * t1041;
t942 = (-t1004 - t1010) * t936 * t830;
t941 = (-t1002 - t1008) * t936 * t831;
t940 = (-t1000 - t1006) * t936 * t832;
t937 = 0.1e1 / pkin(2) ^ 2;
t841 = g(1) * t901 + g(2) * t904;
t840 = g(1) * t900 + g(2) * t903;
t839 = g(1) * t899 + g(2) * t902;
t826 = g(3) * t931 + t841 * t925;
t825 = g(3) * t929 + t840 * t923;
t824 = g(3) * t927 + t839 * t921;
t823 = -g(3) * t925 + t841 * t931;
t822 = -g(3) * t923 + t840 * t929;
t821 = -g(3) * t921 + t839 * t927;
t814 = -t1024 * t904 + t1025 * t901;
t813 = -t1024 * t903 + t1025 * t900;
t812 = -t1024 * t902 + t1025 * t899;
t781 = -t993 * t1013 + t952;
t780 = -t993 * t1015 + t953;
t779 = -t993 * t1017 + t954;
t778 = t781 * t1030 + t909 * t961 + t946;
t777 = t780 * t1031 + t907 * t962 + t947;
t776 = t779 * t1032 + t905 * t963 + t948;
t775 = -t924 * t781 * pkin(1) + t930 * t910 * t961 + t943;
t774 = -t922 * t780 * pkin(1) + t928 * t908 * t962 + t944;
t773 = -t920 * t779 * pkin(1) + t926 * t906 * t963 + t945;
t772 = 0.2e1 * t790 + 0.2e1 * t817 + 0.2e1 * t820 + (-t799 * t802 * t1043 - t793 * t811) * t994 + (t940 - 0.2e1 * t1012) * t995 + t986;
t771 = 0.2e1 * t789 + 0.2e1 * t816 + 0.2e1 * t819 + (-t798 * t801 * t1044 - t792 * t810) * t996 + (t941 - 0.2e1 * t1014) * t997 + t987;
t770 = 0.2e1 * t788 + 0.2e1 * t815 + 0.2e1 * t818 + (-t797 * t800 * t1045 - t791 * t809) * t998 + (t942 - 0.2e1 * t1016) * t999 + t988;
t769 = -t790 * t1043 + (t940 - t1012) * t995 + t952 + t986;
t768 = -t789 * t1044 + (t941 - t1014) * t997 + t953 + t987;
t767 = -t788 * t1045 + (t942 - t1016) * t999 + t954 + t988;
t766 = pkin(1) * (t772 * t930 - 0.2e1 * t796 * t980) + t946;
t765 = pkin(1) * (t771 * t928 - 0.2e1 * t795 * t982) + t947;
t764 = pkin(1) * (t770 * t926 - 0.2e1 * t794 * t984) + t948;
t763 = -pkin(1) * (t924 * t772 + (t802 * t937 + t811 * t936) * t930 * t981) + t943;
t762 = -pkin(1) * (t922 * t771 + (t801 * t937 + t810 * t936) * t928 * t983) + t944;
t761 = -pkin(1) * (t920 * t770 + (t800 * t937 + t809 * t936) * t926 * t985) + t945;
t1 = [(-t812 * t902 - t813 * t903 - t814 * t904) * MDP(1) + t1024 * MDP(8) + ((t779 * t973 + t780 * t970 + t781 * t967) * MDP(2) + (t824 * t973 + t825 * t970 + t826 * t967) * MDP(3) + (t821 * t973 + t822 * t970 + t823 * t967) * MDP(4) + (t767 * t973 + t768 * t970 + t769 * t967) * MDP(5) + (t764 * t973 + t765 * t970 + t766 * t967) * MDP(6) + (t761 * t973 + t762 * t970 + t763 * t967) * MDP(7) + ((-t767 * t979 - t768 * t977 - t769 * t975) * MDP(5) + (-t776 * t979 - t777 * t977 - t778 * t975) * MDP(6) + (-t773 * t979 - t774 * t977 - t775 * t975) * MDP(7)) * t936) * t938; (t812 * t899 + t813 * t900 + t814 * t901) * MDP(1) + t1025 * MDP(8) + ((t779 * t972 + t780 * t969 + t781 * t966) * MDP(2) + (t824 * t972 + t825 * t969 + t826 * t966) * MDP(3) + (t821 * t972 + t822 * t969 + t823 * t966) * MDP(4) + (t767 * t972 + t768 * t969 + t769 * t966) * MDP(5) + (t764 * t972 + t765 * t969 + t766 * t966) * MDP(6) + (t761 * t972 + t762 * t969 + t763 * t966) * MDP(7) + ((-t767 * t978 - t768 * t976 - t769 * t974) * MDP(5) + (-t776 * t978 - t777 * t976 - t778 * t974) * MDP(6) + (-t773 * t978 - t774 * t976 - t775 * t974) * MDP(7)) * t936) * t938; (t917 - g(3)) * MDP(8) + ((-t781 * t1013 - t780 * t1015 - t779 * t1017) * MDP(2) + (-t826 * t1013 - t825 * t1015 - t824 * t1017) * MDP(3) + (-t823 * t1013 - t822 * t1015 - t821 * t1017) * MDP(4) + (-t769 * t1013 - t768 * t1015 - t767 * t1017) * MDP(5) + (-t766 * t1013 - t765 * t1015 - t764 * t1017) * MDP(6) + (-t763 * t1013 - t762 * t1015 - t761 * t1017) * MDP(7) + ((t769 * t1018 + t768 * t1019 + t767 * t1020) * MDP(5) + (t778 * t1018 + t777 * t1019 + t776 * t1020) * MDP(6) + (t775 * t1018 + t774 * t1019 + t773 * t1020) * MDP(7)) * t936) * t938;];
tauX  = t1;
