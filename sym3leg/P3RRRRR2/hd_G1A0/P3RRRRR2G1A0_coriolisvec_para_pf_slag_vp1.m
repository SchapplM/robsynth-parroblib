% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRRRR2G1P1A0
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
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:08
% EndTime: 2020-03-09 21:04:10
% DurationCPUTime: 2.70s
% Computational Cost: add. (13212->270), mult. (10713->498), div. (6999->18), fcn. (9165->54), ass. (0->252)
t965 = 2 * pkin(1);
t764 = qJ(1,3) + legFrame(3,3);
t757 = qJ(2,3) + t764;
t751 = qJ(3,3) + t757;
t752 = -qJ(3,3) + t757;
t955 = -2 * pkin(1);
t706 = sin(t764) * t955 + (-sin(t752) - sin(t751)) * pkin(2);
t801 = qJ(2,3) + qJ(3,3);
t767 = sin(t801);
t802 = qJ(2,3) - qJ(3,3);
t768 = sin(t802);
t734 = 0.1e1 / (t767 + t768);
t827 = xDP(2);
t836 = 0.1e1 / pkin(2);
t839 = 1 / pkin(1);
t903 = t836 * t839;
t885 = t827 * t903;
t691 = t706 * t734 * t885;
t709 = cos(t764) * t955 + (-cos(t752) - cos(t751)) * pkin(2);
t828 = xDP(1);
t884 = t828 * t903;
t694 = t709 * t734 * t884;
t826 = xDP(3);
t810 = sin(qJ(3,3));
t811 = sin(qJ(2,3));
t785 = 0.1e1 / t811;
t920 = t785 * t839;
t894 = t810 * t920;
t816 = cos(qJ(3,3));
t790 = 0.1e1 / t816 ^ 2;
t916 = t790 * t836;
t817 = cos(qJ(2,3));
t937 = t817 * pkin(1);
t938 = t816 * pkin(2);
t855 = (t937 + t938) * t894 * t916;
t849 = t826 * t855;
t685 = t691 + t694 - t849;
t745 = sin(t757);
t748 = cos(t757);
t789 = 0.1e1 / t816;
t872 = t789 * t894;
t904 = t828 * t839;
t905 = t827 * t839;
t700 = t826 * t872 + (t745 * t905 + t748 * t904) * t785;
t682 = t685 + t700;
t835 = pkin(2) ^ 2;
t788 = t816 ^ 2;
t840 = t816 * t788;
t673 = t835 * t682 * t840;
t676 = t694 / 0.2e1 + t691 / 0.2e1 - t849 / 0.2e1 + t700;
t791 = 0.1e1 / t840;
t838 = pkin(1) ^ 2;
t906 = t826 * t836;
t869 = t789 * t817 * t906;
t909 = t810 * t811;
t888 = t826 * t909;
t941 = pkin(2) * t788;
t901 = t817 * t941;
t661 = (((-pkin(1) * t682 * t909 + t789 * t826) * t816 + pkin(1) * t869) * t791 * t906 + ((t673 + t676 * t901 * t965 + (-pkin(1) * t789 * t888 + t838 * t700) * t816) * t700 + (t673 + (t682 * t901 - t888) * pkin(1)) * t685) * t916) * t920;
t800 = t826 ^ 2;
t911 = t800 * t836;
t667 = ((-t816 * t700 * t937 - t682 * t941) * t789 * t700 - t682 * t685 * t938 - t791 * t911) * t920;
t863 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t833 = rSges(3,2) ^ 2;
t834 = rSges(3,1) ^ 2;
t902 = (t833 + t834);
t883 = 2 * rSges(3,3) ^ 2 + t902;
t953 = -m(3) / 0.2e1;
t859 = t883 * t953 + t863;
t763 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t830 = 0.2e1 * qJ(3,3);
t779 = sin(t830);
t782 = cos(t830);
t742 = (-t833 + t834) * m(3) - Icges(3,1) + Icges(3,2);
t951 = -t742 / 0.2e1;
t958 = t763 * t779 + t782 * t951;
t846 = t859 + t958;
t952 = m(3) * rSges(3,3);
t760 = m(2) * rSges(2,2) - t952;
t825 = m(2) * rSges(2,1);
t852 = ((-t825 + (-rSges(3,1) * t816 + rSges(3,2) * t810) * m(3)) * t817 + t760 * t811) * pkin(1);
t688 = t852 + t846;
t703 = -t810 * Icges(3,5) - Icges(3,6) * t816 + (rSges(3,1) * t810 + rSges(3,2) * t816) * m(3) * (pkin(1) * t811 + rSges(3,3));
t761 = -rSges(3,2) * t952 + Icges(3,6);
t762 = rSges(3,1) * t952 - Icges(3,5);
t843 = -m(2) * t838 - Icges(1,3) + (2 * t838 + t883) * t953 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t863;
t858 = -t760 * t817 - t811 * t825;
t862 = t682 * t869;
t910 = t800 / pkin(2) ^ 2;
t891 = t810 * t910;
t868 = t791 * t891;
t898 = t682 * t906;
t929 = t685 * t676;
t882 = t910 / 0.2e1;
t932 = (t790 * t882 + t929) * t811;
t954 = -0.2e1 * t742;
t964 = (-t761 * t810 - t762 * t816) * t790 * t910 + (t789 * t810 * t816 * t954 + (0.2e1 * t789 - 0.4e1 * t816) * t763) * t898 + (t858 * t929 + ((-rSges(3,1) * t932 - rSges(3,2) * t862) * t816 + (-rSges(3,1) * t862 + rSges(3,2) * t932) * t810) * m(3)) * t965 + (t843 + 0.2e1 * t852 + t958) * t667 + t688 * t661 - t703 * t868;
t765 = qJ(1,2) + legFrame(2,3);
t758 = qJ(2,2) + t765;
t753 = qJ(3,2) + t758;
t754 = -qJ(3,2) + t758;
t707 = sin(t765) * t955 + (-sin(t754) - sin(t753)) * pkin(2);
t803 = qJ(2,2) + qJ(3,2);
t769 = sin(t803);
t804 = qJ(2,2) - qJ(3,2);
t770 = sin(t804);
t735 = 0.1e1 / (t769 + t770);
t692 = t707 * t735 * t885;
t710 = cos(t765) * t955 + (-cos(t754) - cos(t753)) * pkin(2);
t695 = t710 * t735 * t884;
t812 = sin(qJ(3,2));
t813 = sin(qJ(2,2));
t786 = 0.1e1 / t813;
t919 = t786 * t839;
t893 = t812 * t919;
t818 = cos(qJ(3,2));
t794 = 0.1e1 / t818 ^ 2;
t914 = t794 * t836;
t819 = cos(qJ(2,2));
t935 = t819 * pkin(1);
t936 = t818 * pkin(2);
t854 = (t935 + t936) * t893 * t914;
t848 = t826 * t854;
t686 = t692 + t695 - t848;
t746 = sin(t758);
t749 = cos(t758);
t793 = 0.1e1 / t818;
t871 = t793 * t893;
t701 = t826 * t871 + (t746 * t905 + t749 * t904) * t786;
t683 = t686 + t701;
t792 = t818 ^ 2;
t841 = t818 * t792;
t674 = t835 * t683 * t841;
t677 = t695 / 0.2e1 + t692 / 0.2e1 - t848 / 0.2e1 + t701;
t795 = 0.1e1 / t841;
t867 = t793 * t819 * t906;
t908 = t812 * t813;
t887 = t826 * t908;
t940 = pkin(2) * t792;
t900 = t819 * t940;
t662 = (((-pkin(1) * t683 * t908 + t793 * t826) * t818 + pkin(1) * t867) * t795 * t906 + ((t674 + t677 * t900 * t965 + (-pkin(1) * t793 * t887 + t838 * t701) * t818) * t701 + (t674 + (t683 * t900 - t887) * pkin(1)) * t686) * t914) * t919;
t668 = ((-t818 * t701 * t935 - t683 * t940) * t793 * t701 - t683 * t686 * t936 - t795 * t911) * t919;
t831 = 0.2e1 * qJ(3,2);
t780 = sin(t831);
t783 = cos(t831);
t957 = t763 * t780 + t783 * t951;
t845 = t859 + t957;
t851 = ((-t825 + (-rSges(3,1) * t818 + rSges(3,2) * t812) * m(3)) * t819 + t760 * t813) * pkin(1);
t689 = t851 + t845;
t704 = -t812 * Icges(3,5) - Icges(3,6) * t818 + (rSges(3,1) * t812 + rSges(3,2) * t818) * m(3) * (pkin(1) * t813 + rSges(3,3));
t857 = -t760 * t819 - t813 * t825;
t861 = t683 * t867;
t890 = t812 * t910;
t866 = t795 * t890;
t897 = t683 * t906;
t928 = t686 * t677;
t931 = (t794 * t882 + t928) * t813;
t963 = (-t761 * t812 - t762 * t818) * t794 * t910 + (t793 * t812 * t818 * t954 + (0.2e1 * t793 - 0.4e1 * t818) * t763) * t897 + (t857 * t928 + ((-rSges(3,1) * t931 - rSges(3,2) * t861) * t818 + (-rSges(3,1) * t861 + rSges(3,2) * t931) * t812) * m(3)) * t965 + (t843 + 0.2e1 * t851 + t957) * t668 + t689 * t662 - t704 * t866;
t766 = qJ(1,1) + legFrame(1,3);
t759 = qJ(2,1) + t766;
t755 = qJ(3,1) + t759;
t756 = -qJ(3,1) + t759;
t708 = sin(t766) * t955 + (-sin(t756) - sin(t755)) * pkin(2);
t805 = qJ(2,1) + qJ(3,1);
t771 = sin(t805);
t806 = qJ(2,1) - qJ(3,1);
t772 = sin(t806);
t736 = 0.1e1 / (t771 + t772);
t693 = t708 * t736 * t885;
t711 = cos(t766) * t955 + (-cos(t756) - cos(t755)) * pkin(2);
t696 = t711 * t736 * t884;
t814 = sin(qJ(3,1));
t815 = sin(qJ(2,1));
t787 = 0.1e1 / t815;
t918 = t787 * t839;
t892 = t814 * t918;
t820 = cos(qJ(3,1));
t798 = 0.1e1 / t820 ^ 2;
t912 = t798 * t836;
t821 = cos(qJ(2,1));
t933 = t821 * pkin(1);
t934 = t820 * pkin(2);
t853 = (t933 + t934) * t892 * t912;
t847 = t826 * t853;
t687 = t693 + t696 - t847;
t747 = sin(t759);
t750 = cos(t759);
t797 = 0.1e1 / t820;
t870 = t797 * t892;
t702 = t826 * t870 + (t747 * t905 + t750 * t904) * t787;
t684 = t687 + t702;
t796 = t820 ^ 2;
t842 = t820 * t796;
t675 = t835 * t684 * t842;
t678 = t696 / 0.2e1 + t693 / 0.2e1 - t847 / 0.2e1 + t702;
t799 = 0.1e1 / t842;
t865 = t797 * t821 * t906;
t907 = t814 * t815;
t886 = t826 * t907;
t939 = pkin(2) * t796;
t899 = t821 * t939;
t663 = (((-pkin(1) * t684 * t907 + t797 * t826) * t820 + pkin(1) * t865) * t799 * t906 + ((t675 + t678 * t899 * t965 + (-pkin(1) * t797 * t886 + t838 * t702) * t820) * t702 + (t675 + (t684 * t899 - t886) * pkin(1)) * t687) * t912) * t918;
t669 = ((-t820 * t702 * t933 - t684 * t939) * t797 * t702 - t684 * t687 * t934 - t799 * t911) * t918;
t832 = 0.2e1 * qJ(3,1);
t781 = sin(t832);
t784 = cos(t832);
t956 = t763 * t781 + t784 * t951;
t844 = t859 + t956;
t850 = ((-t825 + (-rSges(3,1) * t820 + rSges(3,2) * t814) * m(3)) * t821 + t760 * t815) * pkin(1);
t690 = t850 + t844;
t705 = -t814 * Icges(3,5) - Icges(3,6) * t820 + (rSges(3,1) * t814 + rSges(3,2) * t820) * m(3) * (pkin(1) * t815 + rSges(3,3));
t856 = -t760 * t821 - t815 * t825;
t860 = t684 * t865;
t889 = t814 * t910;
t864 = t799 * t889;
t896 = t684 * t906;
t927 = t687 * t678;
t930 = (t798 * t882 + t927) * t815;
t962 = (-t761 * t814 - t762 * t820) * t798 * t910 + (t797 * t814 * t820 * t954 + (0.2e1 * t797 - 0.4e1 * t820) * t763) * t896 + (t856 * t927 + ((-rSges(3,1) * t930 - rSges(3,2) * t860) * t820 + (-rSges(3,1) * t860 + rSges(3,2) * t930) * t814) * m(3)) * t965 + (t843 + 0.2e1 * t850 + t956) * t669 + t690 * t663 - t705 * t864;
t712 = -t761 * t816 + t762 * t810;
t774 = cos(t802);
t895 = t762 * t910;
t923 = t763 * t782;
t926 = t742 * t779;
t944 = pkin(1) * t700 ^ 2;
t947 = cos(t801) / 0.2e1;
t950 = t767 / 0.2e1;
t961 = t846 * t661 + t688 * t667 - t712 * t868 - t790 * t761 * t891 + (-t895 + (-0.2e1 * t923 - t926) * t898) * t789 + (((-t774 / 0.2e1 + t947) * rSges(3,2) + (t768 / 0.2e1 + t950) * rSges(3,1)) * m(3) - t858) * t944;
t713 = -t761 * t818 + t762 * t812;
t776 = cos(t804);
t922 = t763 * t783;
t925 = t742 * t780;
t943 = pkin(1) * t701 ^ 2;
t946 = cos(t803) / 0.2e1;
t949 = t769 / 0.2e1;
t960 = t845 * t662 + t689 * t668 - t713 * t866 - t794 * t761 * t890 + (-t895 + (-0.2e1 * t922 - t925) * t897) * t793 + (((-t776 / 0.2e1 + t946) * rSges(3,2) + (t770 / 0.2e1 + t949) * rSges(3,1)) * m(3) - t857) * t943;
t714 = -t761 * t820 + t762 * t814;
t778 = cos(t806);
t921 = t763 * t784;
t924 = t742 * t781;
t942 = pkin(1) * t702 ^ 2;
t945 = cos(t805) / 0.2e1;
t948 = t771 / 0.2e1;
t959 = t844 * t663 + t690 * t669 - t714 * t864 - t798 * t761 * t889 + (-t895 + (-0.2e1 * t921 - t924) * t896) * t797 + (((-t778 / 0.2e1 + t945) * rSges(3,2) + (t772 / 0.2e1 + t948) * rSges(3,1)) * m(3) - t856) * t942;
t878 = t964 * t785;
t877 = t963 * t786;
t876 = t962 * t787;
t875 = t961 * t734;
t874 = t960 * t735;
t873 = t959 * t736;
t744 = -t902 * m(3) - Icges(3,3);
t1 = [(t750 * t876 + t749 * t877 + t748 * t878 + (t709 * t875 + t710 * t874 + t711 * t873) * t836) * t839; (t747 * t876 + t746 * t877 + t745 * t878 + (t706 * t875 + t707 * t874 + t708 * t873) * t836) * t839; -t959 * t853 - t960 * t854 - t961 * t855 + t962 * t870 + t963 * t871 + t964 * t872 + ((t712 * t661 + t703 * t667 - t744 * t868 + (t926 / 0.2e1 + t923) * t682 ^ 2 + ((t774 / 0.2e1 + t947) * rSges(3,2) + (-t768 / 0.2e1 + t950) * rSges(3,1)) * m(3) * t944) * t789 + (t713 * t662 + t704 * t668 - t744 * t866 + (t925 / 0.2e1 + t922) * t683 ^ 2 + ((t776 / 0.2e1 + t946) * rSges(3,2) + (-t770 / 0.2e1 + t949) * rSges(3,1)) * m(3) * t943) * t793 + (t714 * t663 + t705 * t669 - t744 * t864 + (t924 / 0.2e1 + t921) * t684 ^ 2 + ((t778 / 0.2e1 + t945) * rSges(3,2) + (-t772 / 0.2e1 + t948) * rSges(3,1)) * m(3) * t942) * t797) * t836;];
taucX  = t1;
