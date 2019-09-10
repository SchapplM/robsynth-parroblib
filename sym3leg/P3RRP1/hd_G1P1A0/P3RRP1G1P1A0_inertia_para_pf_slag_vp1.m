% Calculate inertia matrix for parallel robot
% P3RRP1G1P1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:31
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRP1G1P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:03:33
% EndTime: 2019-05-03 15:03:37
% DurationCPUTime: 4.08s
% Computational Cost: add. (24448->471), mult. (51855->673), div. (432->3), fcn. (21590->20), ass. (0->297)
t899 = (qJ(3,3) ^ 2);
t911 = pkin(1) ^ 2;
t890 = cos(qJ(1,3));
t874 = t890 ^ 2;
t952 = t874 * t911;
t1009 = (2 * t899) + t911 - t952;
t900 = (qJ(3,2) ^ 2);
t891 = cos(qJ(1,2));
t875 = t891 ^ 2;
t951 = t875 * t911;
t1008 = (2 * t900) + t911 - t951;
t901 = (qJ(3,1) ^ 2);
t892 = cos(qJ(1,1));
t876 = t892 ^ 2;
t950 = t876 * t911;
t1007 = (2 * t901) + t911 - t950;
t1006 = 2 * pkin(1);
t1005 = 2 * rSges(3,3);
t884 = legFrame(3,3);
t853 = sin(t884);
t833 = t853 * qJ(3,3);
t830 = pkin(2) * t833;
t1004 = -0.2e1 * t830;
t885 = legFrame(2,3);
t854 = sin(t885);
t834 = t854 * qJ(3,2);
t831 = pkin(2) * t834;
t1003 = -0.2e1 * t831;
t886 = legFrame(1,3);
t855 = sin(t886);
t835 = t855 * qJ(3,1);
t832 = pkin(2) * t835;
t1002 = -0.2e1 * t832;
t887 = sin(qJ(1,3));
t1001 = -0.2e1 * t887;
t888 = sin(qJ(1,2));
t1000 = -0.2e1 * t888;
t889 = sin(qJ(1,1));
t999 = -0.2e1 * t889;
t998 = -0.2e1 * t890;
t997 = -0.2e1 * t891;
t996 = -0.2e1 * t892;
t893 = pkin(2) + rSges(3,1);
t877 = qJ(1,3) + qJ(2,3);
t850 = cos(t877);
t841 = t850 ^ 2;
t910 = pkin(2) ^ 2;
t859 = -t899 + t910;
t844 = 0.1e1 + t859;
t847 = sin(t877);
t880 = 0.2e1 + t911;
t949 = t887 * t890;
t986 = pkin(2) * qJ(3,3);
t912 = -t844 * t874 + 0.2e1 * t949 * t986;
t946 = -0.1e1 / 0.2e1 - t910 / 0.2e1;
t964 = t847 * qJ(3,3);
t974 = qJ(3,3) * t887;
t746 = 0.1e1 / ((pkin(1) * (t844 * t949 + (0.2e1 * t874 - 0.1e1) * t986) * t847 + t974) * t850 * t1006 + pkin(1) * t964 * t998 - t899 * t880 + (0.2e1 * (t899 / 0.2e1 - t912 + t946) * t841 + t912) * t911);
t995 = m(3) * t746;
t878 = qJ(1,2) + qJ(2,2);
t851 = cos(t878);
t842 = t851 ^ 2;
t860 = -t900 + t910;
t845 = 0.1e1 + t860;
t848 = sin(t878);
t948 = t888 * t891;
t987 = pkin(2) * qJ(3,2);
t913 = -t845 * t875 + 0.2e1 * t948 * t987;
t962 = t848 * qJ(3,2);
t976 = qJ(3,2) * t888;
t747 = 0.1e1 / ((pkin(1) * (t845 * t948 + (0.2e1 * t875 - 0.1e1) * t987) * t848 + t976) * t851 * t1006 + pkin(1) * t962 * t997 - t900 * t880 + (0.2e1 * (t900 / 0.2e1 - t913 + t946) * t842 + t913) * t911);
t994 = m(3) * t747;
t879 = qJ(1,1) + qJ(2,1);
t852 = cos(t879);
t843 = t852 ^ 2;
t861 = -t901 + t910;
t846 = 0.1e1 + t861;
t849 = sin(t879);
t947 = t889 * t892;
t988 = pkin(2) * qJ(3,1);
t914 = -t846 * t876 + 0.2e1 * t947 * t988;
t960 = t849 * qJ(3,1);
t978 = qJ(3,1) * t889;
t748 = 0.1e1 / ((pkin(1) * (t846 * t947 + (0.2e1 * t876 - 0.1e1) * t988) * t849 + t978) * t852 * t1006 + pkin(1) * t960 * t996 - t901 * t880 + (0.2e1 * (t901 / 0.2e1 - t914 + t946) * t843 + t914) * t911);
t993 = m(3) * t748;
t782 = (t847 * t887 + t850 * t890) * pkin(1) + t893;
t992 = m(3) * t782;
t783 = (t848 * t888 + t851 * t891) * pkin(1) + t893;
t991 = m(3) * t783;
t784 = (t849 * t889 + t852 * t892) * pkin(1) + t893;
t990 = m(3) * t784;
t989 = m(3) * t893;
t856 = cos(t884);
t985 = pkin(2) * t856;
t857 = cos(t885);
t984 = pkin(2) * t857;
t858 = cos(t886);
t983 = pkin(2) * t858;
t982 = pkin(2) * t890;
t981 = pkin(2) * t891;
t980 = pkin(2) * t892;
t836 = t853 * pkin(2);
t837 = t854 * pkin(2);
t838 = t855 * pkin(2);
t979 = qJ(3,1) * t858;
t977 = qJ(3,2) * t857;
t975 = qJ(3,3) * t856;
t802 = t833 + t985;
t973 = t802 * t887;
t806 = t834 + t984;
t972 = t806 * t888;
t810 = t835 + t983;
t971 = t810 * t889;
t970 = t833 * t850;
t969 = t834 * t851;
t968 = t835 * t852;
t967 = t844 * t856;
t966 = t845 * t857;
t965 = t846 * t858;
t963 = t847 * t850;
t961 = t848 * t851;
t959 = t849 * t852;
t958 = t853 * t847;
t957 = t854 * t848;
t956 = t855 * t849;
t955 = t856 * t850;
t954 = t857 * t851;
t953 = t858 * t852;
t882 = t910 / 0.2e1;
t945 = 0.1e1 / 0.2e1 + t882;
t944 = pkin(2) * t979;
t943 = pkin(2) * t977;
t942 = pkin(2) * t975;
t941 = t856 * t964;
t940 = t857 * t962;
t939 = t858 * t960;
t938 = t911 * t949;
t937 = t911 * t948;
t936 = t911 * t947;
t935 = Icges(3,2) + Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t934 = -(2 * t899) - t952;
t933 = -(2 * t900) - t951;
t932 = -(2 * t901) - t950;
t931 = t856 * t938;
t930 = t857 * t937;
t929 = t858 * t936;
t928 = (rSges(3,3) ^ 2) + t910 + (0.2e1 * pkin(2) + rSges(3,1)) * rSges(3,1);
t926 = qJ(3,3) * t1005 + t899 + t928;
t770 = t926 * m(3) + t935;
t925 = qJ(3,2) * t1005 + t900 + t928;
t771 = t925 * m(3) + t935;
t924 = qJ(3,1) * t1005 + t901 + t928;
t772 = t924 * m(3) + t935;
t927 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + m(2) * t911 + t935;
t803 = t833 - t985;
t804 = t836 + t975;
t923 = -t803 * t874 - t804 * t949;
t922 = -t803 * t949 + t804 * t874;
t807 = t834 - t984;
t808 = t837 + t977;
t921 = -t807 * t875 - t808 * t948;
t920 = -t807 * t948 + t808 * t875;
t811 = t835 - t983;
t812 = t838 + t979;
t919 = -t811 * t876 - t812 * t947;
t918 = -t811 * t947 + t812 * t876;
t894 = m(2) * rSges(2,2);
t820 = (-qJ(3,3) - rSges(3,3)) * m(3) + t894;
t829 = m(2) * rSges(2,1) + t989;
t917 = (-(t820 * t890 - t829 * t887) * t847 + (t820 * t887 + t829 * t890) * t850) * pkin(1);
t821 = (-qJ(3,2) - rSges(3,3)) * m(3) + t894;
t916 = (-(t821 * t891 - t829 * t888) * t848 + (t821 * t888 + t829 * t891) * t851) * pkin(1);
t822 = (-qJ(3,1) - rSges(3,3)) * m(3) + t894;
t915 = (-(t822 * t892 - t829 * t889) * t849 + (t822 * t889 + t829 * t892) * t852) * pkin(1);
t907 = koppelP(1,1);
t906 = koppelP(2,1);
t905 = koppelP(3,1);
t904 = koppelP(1,2);
t903 = koppelP(2,2);
t902 = koppelP(3,2);
t898 = rSges(4,1);
t897 = rSges(4,2);
t896 = xP(3);
t881 = 0.1e1 + t911;
t873 = -t901 / 0.2e1;
t872 = -t900 / 0.2e1;
t871 = -t899 / 0.2e1;
t866 = cos(t896);
t865 = sin(t896);
t828 = -0.2e1 * t944;
t827 = -0.2e1 * t943;
t826 = -0.2e1 * t942;
t825 = 0.2e1 * t832;
t824 = 0.2e1 * t831;
t823 = 0.2e1 * t830;
t819 = -t978 + t980;
t818 = -t976 + t981;
t817 = -t974 + t982;
t816 = t855 * t936;
t815 = t854 * t937;
t814 = t853 * t938;
t813 = t838 - t979;
t809 = t837 - t977;
t805 = t836 - t975;
t801 = -t865 * t904 + t866 * t907;
t800 = -t865 * t903 + t866 * t906;
t799 = -t865 * t902 + t866 * t905;
t798 = -t865 * t907 - t866 * t904;
t797 = -t865 * t906 - t866 * t903;
t796 = -t865 * t905 - t866 * t902;
t792 = m(4) * (-t865 * t897 + t866 * t898);
t791 = m(4) * (-t865 * t898 - t866 * t897);
t790 = t855 * t861 + t828;
t789 = t854 * t860 + t827;
t788 = t853 * t859 + t826;
t787 = t846 * t855 + t828;
t786 = t845 * t854 + t827;
t785 = t844 * t853 + t826;
t781 = t787 * t892;
t780 = t786 * t891;
t779 = t785 * t890;
t769 = t813 * t892 - t971;
t768 = t810 * t892 + t813 * t889;
t767 = t809 * t891 - t972;
t766 = t806 * t891 + t809 * t888;
t765 = t805 * t890 - t973;
t764 = t802 * t890 + t805 * t887;
t763 = (t858 * t861 + t825) * t892 + t790 * t889;
t762 = (t857 * t860 + t824) * t891 + t789 * t888;
t761 = (t856 * t859 + t823) * t890 + t788 * t887;
t760 = (t825 + t965) * t892 + t889 * t787;
t759 = (t824 + t966) * t891 + t888 * t786;
t758 = (t823 + t967) * t890 + t887 * t785;
t757 = t790 * t892 + ((t882 + t873) * t858 + t832) * t999;
t756 = t789 * t891 + ((t882 + t872) * t857 + t831) * t1000;
t755 = t788 * t890 + ((t882 + t871) * t856 + t830) * t1001;
t754 = t915 + t772;
t753 = t916 + t771;
t752 = t917 + t770;
t751 = (t911 + t924) * m(3) + 0.2e1 * t915 + t927;
t750 = (t911 + t925) * m(3) + 0.2e1 * t916 + t927;
t749 = (t911 + t926) * m(3) + 0.2e1 * t917 + t927;
t745 = (-t953 + t956) * qJ(3,1) + (-t781 * t843 - t760 * t959 + (t855 * t910 + t855 - t944) * t892 + (-(t1002 - t965) * t843 - qJ(3,1) * t813) * t889) * pkin(1);
t744 = (-t954 + t957) * qJ(3,2) + (-t780 * t842 - t759 * t961 + (t854 * t910 + t854 - t943) * t891 + (-(t1003 - t966) * t842 - qJ(3,2) * t809) * t888) * pkin(1);
t743 = (-t955 + t958) * qJ(3,3) + (-t779 * t841 - t758 * t963 + (t853 * t910 + t853 - t942) * t890 + (-(t1004 - t967) * t841 - qJ(3,3) * t805) * t887) * pkin(1);
t742 = -t939 - t968 + (t760 * t843 - (t781 + ((t873 + t945) * t858 + t832) * t999) * t959 - (t858 * t910 + t832 + t858) * t892 + qJ(3,1) * t971) * pkin(1);
t741 = -t940 - t969 + (t759 * t842 - (t780 + ((t872 + t945) * t857 + t831) * t1000) * t961 - (t857 * t910 + t831 + t857) * t891 + qJ(3,2) * t972) * pkin(1);
t740 = -t941 - t970 + (t758 * t841 - (t779 + ((t871 + t945) * t856 + t830) * t1001) * t963 - (t856 * t910 + t830 + t856) * t890 + qJ(3,3) * t973) * pkin(1);
t739 = (-t1007 * t855 + t828 + t929) * t852 + (t858 * t932 + t816 + t825) * t849 + (-t769 * t843 - t768 * t959 + (t838 - 0.2e1 * t979) * t892 + t889 * t835) * pkin(1);
t738 = (-t1008 * t854 + t827 + t930) * t851 + (t857 * t933 + t815 + t824) * t848 + (-t767 * t842 - t766 * t961 + (t837 - 0.2e1 * t977) * t891 + t888 * t834) * pkin(1);
t737 = (-t1009 * t853 + t826 + t931) * t850 + (t856 * t934 + t814 + t823) * t847 + (-t765 * t841 - t764 * t963 + (t836 - 0.2e1 * t975) * t890 + t887 * t833) * pkin(1);
t736 = (t1007 * t858 + t1002 + t816) * t852 + (t855 * t932 + t828 - t929) * t849 + (-t769 * t959 + t768 * t843 + t835 * t996 + (-t978 - t980) * t858) * pkin(1);
t735 = (t1008 * t857 + t1003 + t815) * t851 + (t854 * t933 + t827 - t930) * t848 + (-t767 * t961 + t766 * t842 + t834 * t997 + (-t976 - t981) * t857) * pkin(1);
t734 = (t1009 * t856 + t1004 + t814) * t850 + (t853 * t934 + t826 - t931) * t847 + (-t765 * t963 + t764 * t841 + t833 * t998 + (-t974 - t982) * t856) * pkin(1);
t733 = -t881 * t939 - t968 + (t757 * t959 - t763 * t843 + t810 * t819) * pkin(1) + ((t919 - t983) * t852 + t918 * t849) * t911;
t732 = -t881 * t940 - t969 + (t756 * t961 - t762 * t842 + t806 * t818) * pkin(1) + ((t921 - t984) * t851 + t920 * t848) * t911;
t731 = -t881 * t941 - t970 + (t755 * t963 - t761 * t841 + t802 * t817) * pkin(1) + ((t923 - t985) * t850 + t922 * t847) * t911;
t730 = (t881 * t956 - t953) * qJ(3,1) + (t757 * t843 + t763 * t959 - t813 * t819) * pkin(1) + ((-t918 + t838) * t852 + t919 * t849) * t911;
t729 = (t881 * t957 - t954) * qJ(3,2) + (t756 * t842 + t762 * t961 - t809 * t818) * pkin(1) + ((-t920 + t837) * t851 + t921 * t848) * t911;
t728 = (t881 * t958 - t955) * qJ(3,3) + (t755 * t841 + t761 * t963 - t805 * t817) * pkin(1) + ((-t922 + t836) * t850 + t923 * t847) * t911;
t727 = (t742 * t801 + t745 * t798) * t748;
t726 = (t741 * t800 + t744 * t797) * t747;
t725 = (t740 * t799 + t743 * t796) * t746;
t724 = (t736 * t801 + t739 * t798) * t748;
t723 = (t735 * t800 + t738 * t797) * t747;
t722 = (t734 * t799 + t737 * t796) * t746;
t721 = (t730 * t798 + t733 * t801) * t748;
t720 = (t729 * t797 + t732 * t800) * t747;
t719 = (t728 * t796 + t731 * t799) * t746;
t718 = (-t730 * t893 - t745 * t784 + t739) * t993;
t717 = (-t729 * t893 - t744 * t783 + t738) * t994;
t716 = (-t728 * t893 - t743 * t782 + t737) * t995;
t715 = (-t733 * t893 - t742 * t784 + t736) * t993;
t714 = (-t732 * t893 - t741 * t783 + t735) * t994;
t713 = (-t731 * t893 - t740 * t782 + t734) * t995;
t712 = (t730 * t772 - t739 * t989 + t745 * t754) * t748;
t711 = (t729 * t771 - t738 * t989 + t744 * t753) * t747;
t710 = (t728 * t770 - t737 * t989 + t743 * t752) * t746;
t709 = (t733 * t772 - t736 * t989 + t742 * t754) * t748;
t708 = (t732 * t771 - t735 * t989 + t741 * t753) * t747;
t707 = (t731 * t770 - t734 * t989 + t740 * t752) * t746;
t706 = (t730 * t754 - t739 * t990 + t745 * t751) * t748;
t705 = (t729 * t753 - t738 * t991 + t744 * t750) * t747;
t704 = (t728 * t752 - t737 * t992 + t743 * t749) * t746;
t703 = (t733 * t754 - t736 * t990 + t742 * t751) * t748;
t702 = (t732 * t753 - t735 * t991 + t741 * t750) * t747;
t701 = (t731 * t752 - t734 * t992 + t740 * t749) * t746;
t700 = (-t721 * t893 - t727 * t784 + t724) * m(3);
t699 = (-t720 * t893 - t726 * t783 + t723) * m(3);
t698 = (-t719 * t893 - t725 * t782 + t722) * m(3);
t697 = t721 * t772 - t724 * t989 + t727 * t754;
t696 = t720 * t771 - t723 * t989 + t726 * t753;
t695 = t719 * t770 - t722 * t989 + t725 * t752;
t694 = t721 * t754 - t724 * t990 + t727 * t751;
t693 = t720 * t753 - t723 * t991 + t726 * t750;
t692 = t719 * t752 - t722 * t992 + t725 * t749;
t1 = [m(4) + (t706 * t745 + t712 * t730 + t718 * t739) * t748 + (t705 * t744 + t711 * t729 + t717 * t738) * t747 + (t704 * t743 + t710 * t728 + t716 * t737) * t746, (t706 * t742 + t712 * t733 + t718 * t736) * t748 + (t705 * t741 + t711 * t732 + t717 * t735) * t747 + (t704 * t740 + t710 * t731 + t716 * t734) * t746, t704 * t725 + t705 * t726 + t706 * t727 + t710 * t719 + t711 * t720 + t712 * t721 + t716 * t722 + t717 * t723 + t718 * t724 + t791; (t703 * t745 + t709 * t730 + t715 * t739) * t748 + (t702 * t744 + t708 * t729 + t714 * t738) * t747 + (t701 * t743 + t707 * t728 + t713 * t737) * t746, m(4) + (t703 * t742 + t709 * t733 + t715 * t736) * t748 + (t702 * t741 + t708 * t732 + t714 * t735) * t747 + (t701 * t740 + t707 * t731 + t713 * t734) * t746, t701 * t725 + t702 * t726 + t703 * t727 + t707 * t719 + t708 * t720 + t709 * t721 + t713 * t722 + t714 * t723 + t715 * t724 + t792; t791 + (t694 * t745 + t697 * t730 + t700 * t739) * t748 + (t693 * t744 + t696 * t729 + t699 * t738) * t747 + (t692 * t743 + t695 * t728 + t698 * t737) * t746, t792 + (t694 * t742 + t697 * t733 + t700 * t736) * t748 + (t693 * t741 + t696 * t732 + t699 * t735) * t747 + (t692 * t740 + t695 * t731 + t698 * t734) * t746, t694 * t727 + t697 * t721 + t700 * t724 + t693 * t726 + t696 * t720 + t699 * t723 + t692 * t725 + t695 * t719 + t698 * t722 + Icges(4,3) + m(4) * (t897 ^ 2 + t898 ^ 2);];
MX  = t1;
