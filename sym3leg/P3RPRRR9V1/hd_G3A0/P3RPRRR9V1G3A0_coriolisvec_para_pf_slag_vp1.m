% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR9V1G3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:56:23
% EndTime: 2020-08-06 18:56:26
% DurationCPUTime: 3.38s
% Computational Cost: add. (18744->360), mult. (24129->591), div. (3921->10), fcn. (20250->40), ass. (0->285)
t1046 = 2 * pkin(1);
t874 = cos(pkin(7));
t1030 = 0.2e1 * t874 ^ 2;
t884 = cos(qJ(3,3));
t869 = t884 ^ 2;
t885 = cos(qJ(1,3));
t1007 = pkin(5) + qJ(2,3);
t859 = pkin(6) + t1007;
t879 = sin(qJ(1,3));
t966 = pkin(1) * t885 + t879 * t859;
t873 = sin(pkin(7));
t878 = sin(qJ(3,3));
t981 = t873 * t878;
t750 = t966 * t981 + (t869 - 0.1e1) * t885 * pkin(3);
t1012 = t884 * pkin(2);
t1015 = t869 * pkin(3);
t765 = pkin(1) * t878 + (-pkin(3) + t1012 + 0.2e1 * t1015) * t873;
t906 = -pkin(3) / 0.2e1;
t785 = t1015 + t1012 / 0.2e1 + t906;
t875 = legFrame(3,2);
t844 = sin(t875);
t847 = cos(t875);
t1032 = 0.2e1 * pkin(3);
t950 = t885 * t981;
t923 = pkin(2) * t950 + (t950 * t1032 - t966) * t884;
t840 = pkin(1) * t873;
t978 = t884 * (-t878 * pkin(3) + t840);
t987 = t844 * t885;
t907 = pkin(2) / 0.2e1;
t994 = (t884 * pkin(3) + t907) * t878;
t729 = (-t785 * t987 + t847 * t994) * t1030 + (t847 * t765 + t923 * t844) * t874 + t750 * t844 + t847 * t978;
t896 = xDP(2);
t778 = 0.1e1 / (t874 * t884 - t981);
t841 = 0.1e1 / t859;
t997 = t778 * t841;
t723 = t729 * t896 * t997;
t984 = t847 * t885;
t730 = (t785 * t984 + t844 * t994) * t1030 + (t844 * t765 - t923 * t847) * t874 - t750 * t847 + t844 * t978;
t897 = xDP(1);
t724 = t730 * t897 * t997;
t866 = pkin(7) + qJ(3,3);
t836 = cos(t866);
t1021 = pkin(3) * t836;
t839 = t874 * pkin(2);
t818 = t839 + pkin(1);
t782 = t818 + t1021;
t762 = -t782 * t879 + t859 * t885;
t895 = xDP(3);
t747 = t762 * t841 * t895;
t720 = t724 + t723 + t747;
t1045 = 0.2e1 * t720;
t886 = cos(qJ(3,2));
t870 = t886 ^ 2;
t887 = cos(qJ(1,2));
t1008 = pkin(5) + qJ(2,2);
t860 = pkin(6) + t1008;
t881 = sin(qJ(1,2));
t965 = pkin(1) * t887 + t881 * t860;
t880 = sin(qJ(3,2));
t980 = t873 * t880;
t751 = t965 * t980 + (t870 - 0.1e1) * t887 * pkin(3);
t1011 = t886 * pkin(2);
t1014 = t870 * pkin(3);
t766 = pkin(1) * t880 + (-pkin(3) + t1011 + 0.2e1 * t1014) * t873;
t786 = t1014 + t1011 / 0.2e1 + t906;
t876 = legFrame(2,2);
t845 = sin(t876);
t848 = cos(t876);
t949 = t887 * t980;
t922 = pkin(2) * t949 + (t949 * t1032 - t965) * t886;
t977 = t886 * (-t880 * pkin(3) + t840);
t986 = t845 * t887;
t993 = (t886 * pkin(3) + t907) * t880;
t731 = (-t786 * t986 + t848 * t993) * t1030 + (t848 * t766 + t922 * t845) * t874 + t751 * t845 + t848 * t977;
t779 = 0.1e1 / (t874 * t886 - t980);
t842 = 0.1e1 / t860;
t996 = t779 * t842;
t725 = t731 * t896 * t996;
t983 = t848 * t887;
t732 = (t786 * t983 + t845 * t993) * t1030 + (t845 * t766 - t922 * t848) * t874 - t751 * t848 + t845 * t977;
t726 = t732 * t897 * t996;
t867 = pkin(7) + qJ(3,2);
t837 = cos(t867);
t1020 = pkin(3) * t837;
t783 = t818 + t1020;
t763 = -t783 * t881 + t860 * t887;
t748 = t763 * t842 * t895;
t721 = t726 + t725 + t748;
t1044 = 0.2e1 * t721;
t888 = cos(qJ(3,1));
t871 = t888 ^ 2;
t889 = cos(qJ(1,1));
t1009 = pkin(5) + qJ(2,1);
t861 = pkin(6) + t1009;
t883 = sin(qJ(1,1));
t964 = pkin(1) * t889 + t883 * t861;
t882 = sin(qJ(3,1));
t979 = t873 * t882;
t752 = t964 * t979 + (t871 - 0.1e1) * t889 * pkin(3);
t1010 = t888 * pkin(2);
t1013 = t871 * pkin(3);
t767 = pkin(1) * t882 + (-pkin(3) + t1010 + 0.2e1 * t1013) * t873;
t787 = t1013 + t1010 / 0.2e1 + t906;
t877 = legFrame(1,2);
t846 = sin(t877);
t849 = cos(t877);
t948 = t889 * t979;
t921 = pkin(2) * t948 + (t948 * t1032 - t964) * t888;
t976 = t888 * (-t882 * pkin(3) + t840);
t985 = t846 * t889;
t992 = (t888 * pkin(3) + t907) * t882;
t733 = (-t787 * t985 + t849 * t992) * t1030 + (t849 * t767 + t921 * t846) * t874 + t752 * t846 + t849 * t976;
t780 = 0.1e1 / (t874 * t888 - t979);
t843 = 0.1e1 / t861;
t995 = t780 * t843;
t727 = t733 * t896 * t995;
t982 = t849 * t889;
t734 = (t787 * t982 + t846 * t992) * t1030 + (t846 * t767 - t921 * t849) * t874 - t752 * t849 + t846 * t976;
t728 = t734 * t897 * t995;
t868 = pkin(7) + qJ(3,1);
t838 = cos(t868);
t1019 = pkin(3) * t838;
t784 = t818 + t1019;
t764 = -t784 * t883 + t861 * t889;
t749 = t764 * t843 * t895;
t722 = t728 + t727 + t749;
t1043 = 0.2e1 * t722;
t830 = sin(t866);
t768 = t847 * t830 - t836 * t987;
t769 = t844 * t830 + t836 * t984;
t824 = 0.1e1 / t836;
t741 = (-t879 * t895 + (t768 * t896 + t769 * t897) * t824) * t841;
t1042 = t741 / 0.2e1;
t831 = sin(t867);
t770 = t848 * t831 - t837 * t986;
t771 = t845 * t831 + t837 * t983;
t825 = 0.1e1 / t837;
t742 = (-t881 * t895 + (t770 * t896 + t771 * t897) * t825) * t842;
t1041 = t742 / 0.2e1;
t832 = sin(t868);
t772 = t849 * t832 - t838 * t985;
t773 = t846 * t832 + t838 * t982;
t826 = 0.1e1 / t838;
t743 = (-t883 * t895 + (t772 * t896 + t773 * t897) * t826) * t843;
t1040 = t743 / 0.2e1;
t1039 = -0.2e1 * t839;
t1038 = t1046 / 0.2e1;
t963 = 2 * m(3);
t1037 = (rSges(3,1) * t832 + rSges(3,2) * t838) * t963;
t1036 = (rSges(3,1) * t831 + rSges(3,2) * t837) * t963;
t1035 = (rSges(3,1) * t830 + rSges(3,2) * t836) * t963;
t1034 = -2 * pkin(1);
t1031 = 4 * rSges(2,3);
t1029 = -4 * pkin(5) - 4 * pkin(6);
t1028 = -m(2) / 0.2e1;
t1027 = pkin(2) * m(3);
t1026 = -rSges(3,1) / 0.2e1;
t1025 = -rSges(3,2) / 0.2e1;
t1024 = m(3) * rSges(3,2);
t1023 = rSges(2,2) * m(2);
t911 = rSges(3,2) ^ 2;
t913 = rSges(3,1) ^ 2;
t796 = (-t911 + t913) * m(3) - Icges(3,1) + Icges(3,2);
t1022 = -t796 / 0.2e1;
t855 = (rSges(3,3) + t1007);
t1018 = t855 * m(3);
t856 = (rSges(3,3) + t1008);
t1017 = t856 * m(3);
t857 = (rSges(3,3) + t1009);
t1016 = t857 * m(3);
t819 = 0.2e1 * t866;
t805 = sin(t819);
t1006 = t741 * t805;
t820 = 0.2e1 * t867;
t806 = sin(t820);
t1005 = t742 * t806;
t821 = 0.2e1 * t868;
t807 = sin(t821);
t1004 = t743 * t807;
t916 = 0.1e1 / pkin(3);
t759 = (t844 * t897 + t847 * t896) * t916 * t824;
t1003 = t759 ^ 2 * t824;
t760 = (t845 * t897 + t848 * t896) * t916 * t825;
t1002 = t760 ^ 2 * t825;
t761 = (t846 * t897 + t849 * t896) * t916 * t826;
t1001 = t761 ^ 2 * t826;
t1000 = t759 * t830;
t999 = t760 * t831;
t998 = t761 * t832;
t991 = (m(2) * rSges(2,1) + t1027) * t874;
t808 = cos(t819);
t823 = rSges(3,1) * t1024 - Icges(3,4);
t990 = t823 * t808;
t809 = cos(t820);
t989 = t823 * t809;
t810 = cos(t821);
t988 = t823 * t810;
t738 = pkin(1) * t741;
t714 = t738 - t724 / 0.2e1 - t723 / 0.2e1 - t747 / 0.2e1;
t905 = 0.2e1 * pkin(7);
t863 = t905 + qJ(3,3);
t833 = cos(t863);
t908 = (qJ(2,3) ^ 2);
t915 = pkin(3) ^ 2;
t858 = cos(t905);
t917 = pkin(2) ^ 2;
t918 = pkin(1) ^ 2;
t920 = -(2 * pkin(6) ^ 2) - t917 * t858 - t915 - t917 - (2 * t918) + ((-4 * pkin(6) - 2 * pkin(5)) * pkin(5));
t699 = ((t714 * t1039 + ((qJ(2,3) * t1029) - t915 * t808 - (2 * t908) + t920) * t1042 + (t782 + t1038) * t720) * t741 + ((t859 * t1000 - 0.2e1 * t714 * t836 + (-pkin(2) * t833 - t1012) * t741) * t741 - (-t859 * t1006 / 0.2e1 + t759 * t782) * t824 * t759) * pkin(3)) * t841;
t711 = (-pkin(3) * t1003 + (-t738 + t1045 + (-t839 - t1021) * t741) * t741) * t841;
t899 = m(2) + m(3);
t956 = t873 * t1023;
t924 = pkin(1) * t899 - t956 + t991;
t937 = rSges(3,1) * t836 - rSges(3,2) * t830;
t744 = t937 * m(3) + t924;
t798 = -rSges(3,2) * t1018 + Icges(3,6);
t801 = rSges(3,1) * t1018 - Icges(3,5);
t753 = -t798 * t836 + t801 * t830;
t789 = t1018 + m(2) * (rSges(2,3) + qJ(2,3));
t827 = sin(t863);
t912 = rSges(2,2) ^ 2;
t914 = rSges(2,1) ^ 2;
t919 = -(m(3) * t917 + (-t912 + t914) * m(2) - Icges(2,1) + Icges(2,2)) * t858 / 0.2e1 - (-rSges(2,1) * t1023 + Icges(2,4)) * sin(t905) + t991 * t1034 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t956 * t1046 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t903 = 2 * t918;
t925 = -t903 / 0.2e1 - t911 / 0.2e1 - t913 / 0.2e1 - t917 / 0.2e1;
t938 = (2 * rSges(2,3) ^ 2) + t903 + t912 + t914;
t953 = t830 * t1003;
t954 = t1024 * t1034;
t955 = m(3) * rSges(3,1) * t1046;
t960 = t759 * t1027;
t969 = t833 + t884;
t972 = t827 + t878;
t975 = t744 * t699 - t753 * t953 + (t808 * t1022 + t823 * t805 + ((qJ(2,3) * t1031) + (2 * t908) + t938) * t1028 + t919 + (-(t855 ^ 2) - t937 * t1046 + (-t969 * rSges(3,1) + t972 * rSges(3,2)) * pkin(2) + t925) * m(3)) * t711 + (-t801 * t836 * t759 - t798 * t1000 - t1006 * t796) * t759 + (-t955 * t1000 + t789 * t1045 + 0.2e1 * (t884 * t1025 + t878 * t1026) * t960 + (-rSges(3,1) * t827 - rSges(3,2) * t833) * t960 + (t836 * t954 - 0.2e1 * t990) * t759) * t741;
t739 = pkin(1) * t742;
t715 = t739 - t726 / 0.2e1 - t725 / 0.2e1 - t748 / 0.2e1;
t864 = t905 + qJ(3,2);
t834 = cos(t864);
t909 = (qJ(2,2) ^ 2);
t700 = ((t715 * t1039 + ((qJ(2,2) * t1029) - t915 * t809 - (2 * t909) + t920) * t1041 + (t783 + t1038) * t721) * t742 + ((-0.2e1 * t715 * t837 + t860 * t999 + (-pkin(2) * t834 - t1011) * t742) * t742 - (-t860 * t1005 / 0.2e1 + t760 * t783) * t825 * t760) * pkin(3)) * t842;
t712 = (-pkin(3) * t1002 + (-t739 + t1044 + (-t839 - t1020) * t742) * t742) * t842;
t935 = rSges(3,1) * t837 - rSges(3,2) * t831;
t745 = t935 * m(3) + t924;
t799 = -rSges(3,2) * t1017 + Icges(3,6);
t802 = rSges(3,1) * t1017 - Icges(3,5);
t754 = -t799 * t837 + t802 * t831;
t790 = t1017 + m(2) * (rSges(2,3) + qJ(2,2));
t828 = sin(t864);
t952 = t831 * t1002;
t959 = t760 * t1027;
t968 = t834 + t886;
t971 = t828 + t880;
t974 = t745 * t700 - t754 * t952 + (t809 * t1022 + t823 * t806 + ((qJ(2,2) * t1031) + (2 * t909) + t938) * t1028 + t919 + (-(t856 ^ 2) - t935 * t1046 + (-t968 * rSges(3,1) + t971 * rSges(3,2)) * pkin(2) + t925) * m(3)) * t712 + (-t802 * t837 * t760 - t1005 * t796 - t799 * t999) * t760 + (-t955 * t999 + t790 * t1044 + 0.2e1 * (t886 * t1025 + t880 * t1026) * t959 + (-rSges(3,1) * t828 - rSges(3,2) * t834) * t959 + (t837 * t954 - 0.2e1 * t989) * t760) * t742;
t740 = pkin(1) * t743;
t716 = t740 - t728 / 0.2e1 - t727 / 0.2e1 - t749 / 0.2e1;
t865 = t905 + qJ(3,1);
t835 = cos(t865);
t910 = (qJ(2,1) ^ 2);
t701 = ((t716 * t1039 + ((qJ(2,1) * t1029) - t915 * t810 - (2 * t910) + t920) * t1040 + (t784 + t1038) * t722) * t743 + ((-0.2e1 * t716 * t838 + t861 * t998 + (-pkin(2) * t835 - t1010) * t743) * t743 - (-t861 * t1004 / 0.2e1 + t761 * t784) * t826 * t761) * pkin(3)) * t843;
t713 = (-pkin(3) * t1001 + (-t740 + t1043 + (-t839 - t1019) * t743) * t743) * t843;
t933 = rSges(3,1) * t838 - rSges(3,2) * t832;
t746 = t933 * m(3) + t924;
t800 = -rSges(3,2) * t1016 + Icges(3,6);
t803 = rSges(3,1) * t1016 - Icges(3,5);
t755 = -t800 * t838 + t803 * t832;
t791 = t1016 + m(2) * (rSges(2,3) + qJ(2,1));
t829 = sin(t865);
t951 = t832 * t1001;
t958 = t761 * t1027;
t967 = t835 + t888;
t970 = t829 + t882;
t973 = t746 * t701 - t755 * t951 + (t810 * t1022 + t823 * t807 + ((qJ(2,1) * t1031) + (2 * t910) + t938) * t1028 + t919 + (-(t857 ^ 2) - t933 * t1046 + (-t967 * rSges(3,1) + t970 * rSges(3,2)) * pkin(2) + t925) * m(3)) * t713 + (-t803 * t838 * t761 - t1004 * t796 - t800 * t998) * t761 + (-t955 * t998 + t791 * t1043 + 0.2e1 * (t888 * t1025 + t882 * t1026) * t958 + (-rSges(3,1) * t829 - rSges(3,2) * t835) * t958 + (t838 * t954 - 0.2e1 * t988) * t761) * t743;
t944 = t975 * t824;
t943 = t974 * t825;
t942 = t973 * t826;
t941 = -t899 * t699 + t744 * t711 - (-t759 * t1035 + t789 * t741) * t741;
t940 = -t899 * t700 + t745 * t712 - (-t760 * t1036 + t790 * t742) * t742;
t939 = -t899 * t701 + t746 * t713 - (-t761 * t1037 + t791 * t743) * t743;
t931 = t941 * t778;
t930 = t940 * t779;
t929 = t939 * t780;
t804 = -(t911 + t913) * m(3) - Icges(3,3);
t928 = t824 * (t753 * t711 - t804 * t953 + ((t738 - 0.2e1 * t724 - 0.2e1 * t723 - 0.2e1 * t747) * t1035 + (t796 * t805 + 0.2e1 * t990 + (t972 * rSges(3,1) + t969 * rSges(3,2)) * t1027) * t741) * t1042);
t927 = t825 * (t754 * t712 - t804 * t952 + ((t739 - 0.2e1 * t726 - 0.2e1 * t725 - 0.2e1 * t748) * t1036 + (t796 * t806 + 0.2e1 * t989 + (t971 * rSges(3,1) + t968 * rSges(3,2)) * t1027) * t742) * t1041);
t926 = t826 * (t755 * t713 - t804 * t951 + ((t740 - 0.2e1 * t728 - 0.2e1 * t727 - 0.2e1 * t749) * t1037 + (t796 * t807 + 0.2e1 * t988 + (t970 * rSges(3,1) + t967 * rSges(3,2)) * t1027) * t743) * t1040);
t1 = [(t734 * t929 + t773 * t942) * t843 + (t732 * t930 + t771 * t943) * t842 + (t730 * t931 + t769 * t944) * t841 + (t844 * t928 + t845 * t927 + t846 * t926) * t916; (t733 * t929 + t772 * t942) * t843 + (t731 * t930 + t770 * t943) * t842 + (t729 * t931 + t768 * t944) * t841 + (t847 * t928 + t848 * t927 + t849 * t926) * t916; (t939 * t764 - t973 * t883) * t843 + (t940 * t763 - t974 * t881) * t842 + (t941 * t762 - t975 * t879) * t841;];
taucX  = t1;
