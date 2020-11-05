% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4*4x15]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR1G2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_inertia_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:59
% EndTime: 2020-08-07 10:59:16
% DurationCPUTime: 18.02s
% Computational Cost: add. (3437->652), mult. (11532->1377), div. (4620->22), fcn. (13296->26), ass. (0->559)
t567 = sin(qJ(2,1));
t572 = cos(qJ(3,1));
t548 = 0.1e1 / t572;
t566 = sin(qJ(3,1));
t895 = t548 * t566;
t534 = 0.1e1 / t567;
t573 = cos(qJ(2,1));
t553 = t573 ^ 2;
t913 = t534 * t553;
t660 = (t567 + t913) * t895;
t565 = sin(qJ(2,2));
t570 = cos(qJ(3,2));
t542 = 0.1e1 / t570;
t564 = sin(qJ(3,2));
t901 = t542 * t564;
t529 = 0.1e1 / t565;
t571 = cos(qJ(2,2));
t547 = t571 ^ 2;
t921 = t529 * t547;
t661 = (t565 + t921) * t901;
t563 = sin(qJ(2,3));
t568 = cos(qJ(3,3));
t536 = 0.1e1 / t568;
t562 = sin(qJ(3,3));
t907 = t536 * t562;
t524 = 0.1e1 / t563;
t569 = cos(qJ(2,3));
t541 = t569 ^ 2;
t929 = t524 * t541;
t662 = (t563 + t929) * t907;
t555 = sin(qJ(2,4));
t556 = cos(qJ(3,4));
t515 = 0.1e1 / t556;
t554 = sin(qJ(3,4));
t937 = t515 * t554;
t513 = 0.1e1 / t555;
t557 = cos(qJ(2,4));
t520 = t557 ^ 2;
t943 = t513 * t520;
t663 = (t555 + t943) * t937;
t605 = t572 ^ 2;
t549 = 0.1e1 / t605;
t893 = t549 * t566;
t535 = 0.1e1 / t567 ^ 2;
t909 = t535 * t553;
t664 = (0.1e1 + t909) * t893;
t601 = t570 ^ 2;
t543 = 0.1e1 / t601;
t899 = t543 * t564;
t530 = 0.1e1 / t565 ^ 2;
t917 = t530 * t547;
t665 = (0.1e1 + t917) * t899;
t597 = t568 ^ 2;
t537 = 0.1e1 / t597;
t905 = t537 * t562;
t525 = 0.1e1 / t563 ^ 2;
t925 = t525 * t541;
t666 = (0.1e1 + t925) * t905;
t587 = t556 ^ 2;
t516 = 0.1e1 / t587;
t935 = t516 * t554;
t514 = 0.1e1 / t555 ^ 2;
t939 = t514 * t520;
t667 = (0.1e1 + t939) * t935;
t574 = xP(4);
t508 = sin(t574);
t509 = cos(t574);
t578 = koppelP(1,2);
t582 = koppelP(1,1);
t486 = t508 * t582 + t509 * t578;
t490 = -t508 * t578 + t509 * t582;
t561 = legFrame(1,2);
t503 = sin(t561);
t507 = cos(t561);
t882 = t567 * t572;
t482 = t503 * t566 + t507 * t882;
t915 = t534 * t548;
t858 = t482 * t915;
t946 = t507 * t566;
t479 = t503 * t882 - t946;
t861 = t479 * t915;
t680 = -t486 * t861 + t490 * t858;
t577 = koppelP(2,2);
t581 = koppelP(2,1);
t485 = t508 * t581 + t509 * t577;
t489 = -t508 * t577 + t509 * t581;
t560 = legFrame(2,2);
t502 = sin(t560);
t506 = cos(t560);
t884 = t565 * t570;
t481 = t502 * t564 + t506 * t884;
t923 = t529 * t542;
t859 = t481 * t923;
t948 = t506 * t564;
t478 = t502 * t884 - t948;
t862 = t478 * t923;
t681 = -t485 * t862 + t489 * t859;
t576 = koppelP(3,2);
t580 = koppelP(3,1);
t484 = t508 * t580 + t509 * t576;
t488 = -t508 * t576 + t509 * t580;
t559 = legFrame(3,2);
t501 = sin(t559);
t505 = cos(t559);
t886 = t563 * t568;
t480 = t501 * t562 + t505 * t886;
t931 = t524 * t536;
t860 = t480 * t931;
t950 = t505 * t562;
t477 = t501 * t886 - t950;
t863 = t477 * t931;
t682 = -t484 * t863 + t488 * t860;
t575 = koppelP(4,2);
t579 = koppelP(4,1);
t483 = t508 * t579 + t509 * t575;
t487 = -t508 * t575 + t509 * t579;
t558 = legFrame(4,2);
t500 = sin(t558);
t504 = cos(t558);
t888 = t555 * t556;
t476 = t500 * t554 + t504 * t888;
t945 = t513 * t515;
t864 = t476 * t945;
t952 = t504 * t554;
t475 = t500 * t888 - t952;
t865 = t475 * t945;
t683 = -t483 * t865 + t487 * t864;
t1005 = t515 * t683 + t536 * t682 + t542 * t681 + t548 * t680;
t894 = t548 * t573;
t827 = t534 * t894;
t751 = t566 * t827;
t900 = t542 * t571;
t836 = t529 * t900;
t763 = t564 * t836;
t906 = t536 * t569;
t845 = t524 * t906;
t775 = t562 * t845;
t936 = t515 * t557;
t856 = t513 * t936;
t787 = t554 * t856;
t1004 = t680 * t751 + t681 * t763 + t682 * t775 + t683 * t787;
t1003 = -t680 * t827 - t681 * t836 - t682 * t845 - t683 * t856;
t1002 = t680 * t858 + t681 * t859 + t682 * t860 + t683 * t864;
t1001 = t680 * t861 + t681 * t862 + t682 * t863 + t683 * t865;
t912 = t534 * t573;
t920 = t529 * t571;
t928 = t524 * t569;
t942 = t513 * t557;
t1000 = t680 * t912 + t681 * t920 + t682 * t928 + t683 * t942;
t583 = 1 / pkin(2);
t999 = 2 * t583;
t584 = 1 / (pkin(2) ^ 2);
t998 = 2 * t584;
t444 = t483 * t504 + t487 * t500;
t997 = t683 * t444;
t445 = t484 * t505 + t488 * t501;
t996 = t682 * t445;
t446 = t485 * t506 + t489 * t502;
t995 = t681 * t446;
t443 = t486 * t507 + t490 * t503;
t994 = t680 * t443;
t889 = t554 * t557;
t848 = t516 * t889;
t785 = t513 * t848;
t726 = t444 * t785;
t430 = t583 * t726;
t993 = t430 * t515;
t992 = t430 * t557;
t887 = t562 * t569;
t819 = t537 * t887;
t773 = t524 * t819;
t723 = t445 * t773;
t431 = t583 * t723;
t991 = t431 * t536;
t990 = t431 * t569;
t885 = t564 * t571;
t817 = t543 * t885;
t761 = t529 * t817;
t720 = t446 * t761;
t432 = t583 * t720;
t989 = t432 * t542;
t988 = t432 * t571;
t883 = t566 * t573;
t815 = t549 * t883;
t749 = t534 * t815;
t729 = t443 * t749;
t433 = t583 * t729;
t987 = t433 * t548;
t986 = t433 * t573;
t973 = t444 * t515;
t439 = t583 * t973;
t985 = t439 * t513;
t984 = t439 * t515;
t970 = t445 * t536;
t440 = t583 * t970;
t983 = t440 * t524;
t982 = t440 * t536;
t967 = t446 * t542;
t441 = t583 * t967;
t981 = t441 * t529;
t980 = t441 * t542;
t976 = t443 * t548;
t442 = t583 * t976;
t979 = t442 * t534;
t978 = t442 * t548;
t977 = t443 * t503;
t550 = t548 * t549;
t975 = t443 * t550;
t974 = t444 * t500;
t517 = t515 * t516;
t972 = t444 * t517;
t971 = t445 * t501;
t538 = t536 * t537;
t969 = t445 * t538;
t968 = t446 * t502;
t544 = t542 * t543;
t966 = t446 * t544;
t965 = t475 * t504;
t964 = t476 * t500;
t963 = t477 * t505;
t962 = t478 * t506;
t961 = t479 * t507;
t960 = t480 * t501;
t959 = t481 * t502;
t958 = t482 * t503;
t957 = t500 * t516;
t956 = t501 * t537;
t955 = t502 * t543;
t954 = t503 * t549;
t953 = t504 * t516;
t951 = t505 * t537;
t949 = t506 * t543;
t947 = t507 * t549;
t944 = t513 * t516;
t941 = t514 * t516;
t519 = t557 * t520;
t940 = t514 * t519;
t938 = t514 * t557;
t934 = t517 * t557;
t933 = t519 * t554;
t932 = t520 * t554;
t930 = t524 * t537;
t927 = t525 * t537;
t540 = t569 * t541;
t926 = t525 * t540;
t924 = t525 * t569;
t922 = t529 * t543;
t919 = t530 * t543;
t546 = t571 * t547;
t918 = t530 * t546;
t916 = t530 * t571;
t914 = t534 * t549;
t911 = t535 * t549;
t552 = t573 * t553;
t910 = t535 * t552;
t908 = t535 * t573;
t904 = t538 * t569;
t903 = t540 * t562;
t902 = t541 * t562;
t898 = t544 * t571;
t897 = t546 * t564;
t896 = t547 * t564;
t892 = t550 * t573;
t891 = t552 * t566;
t890 = t553 * t566;
t881 = t683 * t944;
t880 = t682 * t930;
t879 = t681 * t922;
t878 = t680 * t914;
t877 = t430 * t945;
t876 = t431 * t931;
t875 = t432 * t923;
t874 = t433 * t915;
t873 = t534 * t975;
t872 = t535 * t975;
t871 = t513 * t972;
t870 = t514 * t972;
t869 = t524 * t969;
t868 = t525 * t969;
t867 = t529 * t966;
t866 = t530 * t966;
t512 = t554 ^ 2;
t857 = t512 * t941;
t855 = t513 * t935;
t854 = t516 * t942;
t853 = t514 * t936;
t852 = t516 * t940;
t851 = t517 * t939;
t850 = 0.1e1 / t587 ^ 2 * t939;
t849 = t514 * t932;
t847 = t517 * t889;
t523 = t562 ^ 2;
t846 = t523 * t927;
t844 = t524 * t905;
t843 = t537 * t928;
t842 = t525 * t906;
t841 = t537 * t926;
t840 = t538 * t925;
t839 = 0.1e1 / t597 ^ 2 * t925;
t838 = t525 * t902;
t528 = t564 ^ 2;
t837 = t528 * t919;
t835 = t529 * t899;
t834 = t543 * t920;
t833 = t530 * t900;
t832 = t543 * t918;
t831 = t544 * t917;
t830 = 0.1e1 / t601 ^ 2 * t917;
t829 = t530 * t896;
t533 = t566 ^ 2;
t828 = t533 * t911;
t826 = t534 * t893;
t825 = t549 * t912;
t824 = t535 * t894;
t823 = t549 * t910;
t822 = t550 * t909;
t821 = 0.1e1 / t605 ^ 2 * t909;
t820 = t535 * t890;
t818 = t538 * t887;
t816 = t544 * t885;
t814 = t550 * t883;
t808 = t557 + t940;
t806 = t569 + t926;
t804 = t571 + t918;
t802 = t573 + t910;
t801 = t683 * t848;
t800 = t682 * t819;
t799 = t681 * t817;
t798 = t680 * t815;
t797 = t430 * t856;
t796 = t431 * t845;
t795 = t432 * t836;
t794 = t433 * t827;
t510 = t512 ^ 2;
t793 = t510 * t850;
t511 = t554 * t512;
t792 = t511 * t851;
t791 = t511 * t514 * t934;
t790 = t512 * t513 * t934;
t789 = t557 * t857;
t788 = t512 * t850;
t786 = t520 * t855;
t784 = t513 * t847;
t783 = t517 * t849;
t782 = t514 * t847;
t521 = t523 ^ 2;
t781 = t521 * t839;
t522 = t562 * t523;
t780 = t522 * t840;
t779 = t522 * t525 * t904;
t778 = t523 * t524 * t904;
t777 = t569 * t846;
t776 = t523 * t839;
t774 = t541 * t844;
t772 = t524 * t818;
t771 = t538 * t838;
t770 = t525 * t818;
t526 = t528 ^ 2;
t769 = t526 * t830;
t527 = t564 * t528;
t768 = t527 * t831;
t767 = t527 * t530 * t898;
t766 = t528 * t529 * t898;
t765 = t571 * t837;
t764 = t528 * t830;
t762 = t547 * t835;
t760 = t529 * t816;
t759 = t544 * t829;
t758 = t530 * t816;
t531 = t533 ^ 2;
t757 = t531 * t821;
t532 = t566 * t533;
t756 = t532 * t822;
t755 = t532 * t535 * t892;
t754 = t533 * t534 * t892;
t753 = t573 * t828;
t752 = t533 * t821;
t750 = t553 * t826;
t748 = t534 * t814;
t747 = t550 * t820;
t746 = t535 * t814;
t692 = t535 * t583 * t815;
t693 = t530 * t583 * t817;
t694 = t525 * t583 * t819;
t695 = t514 * t583 * t848;
t745 = t475 * t695 + t477 * t694 + t478 * t693 + t479 * t692;
t744 = t476 * t695 + t480 * t694 + t481 * t693 + t482 * t692;
t743 = t683 * t786;
t742 = t682 * t774;
t741 = t681 * t762;
t740 = t680 * t750;
t739 = t430 * t511 * t854;
t738 = t512 * t797;
t737 = t431 * t522 * t843;
t736 = t523 * t796;
t735 = t432 * t527 * t834;
t734 = t528 * t795;
t733 = t433 * t532 * t825;
t732 = t533 * t794;
t731 = t821 * t977;
t730 = t443 * t750;
t728 = t850 * t974;
t727 = t444 * t786;
t725 = t839 * t971;
t724 = t445 * t774;
t722 = t830 * t968;
t721 = t446 * t762;
t719 = t500 * t785;
t718 = t501 * t773;
t717 = t502 * t761;
t716 = t503 * t749;
t715 = t504 * t793;
t714 = t504 * t792;
t713 = t504 * t790;
t712 = t504 * t788;
t711 = t504 * t785;
t710 = t505 * t781;
t709 = t505 * t780;
t708 = t505 * t778;
t707 = t505 * t776;
t706 = t505 * t773;
t705 = t506 * t769;
t704 = t506 * t768;
t703 = t506 * t766;
t702 = t506 * t764;
t701 = t506 * t761;
t700 = t507 * t757;
t699 = t507 * t756;
t698 = t507 * t754;
t697 = t507 * t752;
t696 = t507 * t749;
t691 = t554 * t808;
t690 = t562 * t806;
t689 = t564 * t804;
t688 = t566 * t802;
t687 = t475 * t500 - t476 * t504;
t686 = t477 * t501 - t480 * t505;
t685 = t478 * t502 - t481 * t506;
t684 = t479 * t503 - t482 * t507;
t679 = t512 * t516 * t943 - t555;
t678 = t512 * t851 - t515;
t677 = t512 * t852 - t557;
t676 = t523 * t537 * t929 - t563;
t675 = t523 * t840 - t536;
t674 = t523 * t841 - t569;
t673 = t528 * t543 * t921 - t565;
t672 = t528 * t831 - t542;
t671 = t528 * t832 - t571;
t670 = t533 * t549 * t913 - t567;
t669 = t533 * t822 - t548;
t668 = t533 * t823 - t573;
t659 = t439 * t787 + t430;
t658 = t440 * t775 + t431;
t657 = t441 * t763 + t432;
t656 = t442 * t751 + t433;
t655 = t443 * t664;
t654 = t444 * t667;
t653 = t445 * t666;
t652 = t446 * t665;
t651 = t500 * t667;
t650 = t501 * t666;
t649 = t502 * t665;
t648 = t503 * t664;
t647 = t504 * t667;
t646 = t505 * t666;
t645 = t506 * t665;
t644 = t507 * t664;
t643 = t443 * t669;
t642 = t444 * t678;
t641 = t445 * t675;
t640 = t446 * t672;
t639 = t500 * t679;
t638 = t500 * t678;
t637 = t501 * t676;
t636 = t501 * t675;
t635 = t502 * t673;
t634 = t502 * t672;
t633 = t503 * t670;
t632 = t503 * t669;
t631 = t504 * t679;
t630 = t504 * t678;
t629 = t505 * t676;
t628 = t506 * t673;
t627 = t506 * t672;
t626 = t507 * t670;
t625 = t507 * t669;
t624 = t675 * t505;
t623 = t500 * t663;
t622 = t501 * t662;
t621 = t502 * t661;
t620 = t503 * t660;
t619 = t504 * t663;
t618 = t505 * t662;
t617 = t506 * t661;
t616 = t507 * t660;
t615 = t439 * t512 * t854 + t430 * t937;
t614 = t440 * t523 * t843 + t431 * t907;
t613 = t441 * t528 * t834 + t432 * t901;
t612 = t442 * t533 * t825 + t433 * t895;
t611 = t909 + t917 + t925 + t939;
t409 = t475 * t853 + t477 * t842 + t478 * t833 + t479 * t824;
t410 = t476 * t853 + t480 * t842 + t481 * t833 + t482 * t824;
t499 = t507 ^ 2;
t498 = t506 ^ 2;
t497 = t505 ^ 2;
t496 = t504 ^ 2;
t495 = t503 ^ 2;
t494 = t502 ^ 2;
t493 = t501 ^ 2;
t492 = t500 ^ 2;
t491 = t508 ^ 2 + t509 ^ 2;
t474 = t583 * t626;
t473 = t583 * t633;
t472 = t583 * t628;
t471 = t583 * t635;
t470 = t583 * t629;
t469 = t583 * t637;
t468 = t583 * t631;
t467 = t583 * t639;
t460 = t583 * t616;
t459 = t583 * t620;
t458 = t583 * t617;
t457 = t583 * t621;
t456 = t583 * t618;
t455 = t583 * t622;
t452 = t583 * t619;
t451 = t583 * t623;
t438 = (t504 * t945 + t505 * t931 + t506 * t923 + t507 * t915) * t584;
t437 = (-t500 * t945 - t501 * t931 - t502 * t923 - t503 * t915) * t584;
t436 = (-t500 * t953 - t501 * t951 - t502 * t949 - t503 * t947) * t584;
t435 = (t504 * t855 + t505 * t844 + t506 * t835 + t507 * t826) * t584;
t434 = (-t500 * t855 - t501 * t844 - t502 * t835 - t503 * t826) * t584;
t429 = (-t504 * t782 - t505 * t770 - t506 * t758 - t507 * t746) * t584;
t428 = (-t504 * t791 - t505 * t779 - t506 * t767 - t507 * t755) * t584;
t427 = (t500 * t782 + t501 * t770 + t502 * t758 + t503 * t746) * t584;
t426 = (t500 * t791 + t501 * t779 + t502 * t767 + t503 * t755) * t584;
t425 = (-t504 * t789 - t505 * t777 - t506 * t765 - t507 * t753) * t998;
t424 = (t500 * t789 + t501 * t777 + t502 * t765 + t503 * t753) * t998;
t415 = (-t500 * t712 - t501 * t707 - t502 * t702 - t503 * t697) * t584;
t414 = (-t500 * t715 - t501 * t710 - t502 * t705 - t503 * t700) * t584;
t413 = (t500 * t711 + t501 * t706 + t502 * t701 + t503 * t696) * t998;
t412 = (t500 * t713 + t501 * t708 + t502 * t703 + t503 * t698) * t998;
t411 = (-t500 * t714 - t501 * t709 - t502 * t704 - t503 * t699) * t998;
t408 = t433 * t883 - t442 * t882;
t407 = t432 * t885 - t441 * t884;
t406 = t431 * t887 - t440 * t886;
t405 = -t442 * t566 * t567 - t572 * t986;
t404 = -t441 * t564 * t565 - t570 * t988;
t403 = -t440 * t562 * t563 - t568 * t990;
t402 = t430 * t889 - t439 * t888;
t401 = -t439 * t554 * t555 - t556 * t992;
t400 = ((t503 * t890 + t482) * t914 + (t502 * t896 + t481) * t922 + (t501 * t902 + t480) * t930 + (t500 * t932 + t476) * t944) * t583;
t399 = ((-t507 * t890 + t479) * t914 + (-t506 * t896 + t478) * t922 + (-t505 * t902 + t477) * t930 + (-t504 * t932 + t475) * t944) * t583;
t398 = t475 * t476 * t941 + t477 * t480 * t927 + t478 * t481 * t919 + t479 * t482 * t911;
t397 = ((-t482 * t573 - t503 * t891) * t911 + (-t481 * t571 - t502 * t897) * t919 + (-t480 * t569 - t501 * t903) * t927 + (-t476 * t557 - t500 * t933) * t941) * t583;
t396 = ((-t479 * t573 + t507 * t891) * t911 + (-t478 * t571 + t506 * t897) * t919 + (-t477 * t569 + t505 * t903) * t927 + (-t475 * t557 + t504 * t933) * t941) * t583;
t395 = (t684 * t748 + t685 * t760 + t686 * t772 + t687 * t784) * t583;
t394 = (-t684 * t747 - t685 * t759 - t686 * t771 - t687 * t783) * t583;
t1 = [t475 ^ 2 * t941 + t477 ^ 2 * t927 + t478 ^ 2 * t919 + t479 ^ 2 * t911, (t496 * t788 + t497 * t776 + t498 * t764 + t499 * t752) * t584, (t747 * t961 + t759 * t962 + t771 * t963 + t783 * t965) * t999, (-t748 * t961 - t760 * t962 - t772 * t963 - t784 * t965) * t999, (t496 * t793 + t497 * t781 + t498 * t769 + t499 * t757) * t584, (t496 * t792 + t497 * t780 + t498 * t768 + t499 * t756) * t998, (-t496 * t790 - t497 * t778 - t498 * t766 - t499 * t754) * t998, (-t496 * t785 - t497 * t773 - t498 * t761 - t499 * t749) * t998, (t496 * t516 + t497 * t537 + t498 * t543 + t499 * t549) * t584, t452 * t865 + t456 * t863 + t458 * t862 + t460 * t861 + (t475 * t647 + t477 * t646 + t478 * t645 + t479 * t644) * t583, -t468 * t865 - t470 * t863 - t472 * t862 - t474 * t861 + (-t475 * t630 - t477 * t624 - t478 * t627 - t479 * t625) * t583, 0, 0, 0, t491; t398, t415, t394, t395, t414, t411, t412, t413, t436, -t451 * t865 - t455 * t863 - t457 * t862 - t459 * t861 + (t476 * t647 + t480 * t646 + t481 * t645 + t482 * t644) * t583, t467 * t865 + t469 * t863 + t471 * t862 + t473 * t861 + (-t476 * t630 - t480 * t624 - t481 * t627 - t482 * t625) * t583, 0, 0, 0, 0; t409, t429, t396, t399, t428, t425, t435, t438, 0, ((-t479 * t908 + t802 * t946) * t548 + (-t478 * t916 + t804 * t948) * t542 + (-t477 * t924 + t806 * t950) * t536 + (-t475 * t938 + t808 * t952) * t515) * t583, (-t504 * t677 - t505 * t674 - t506 * t671 - t507 * t668) * t583 + t745, 0, 0, 0, 0; t1001, (-t430 * t711 - t431 * t706 - t432 * t701 - t433 * t696) * t583, -t475 * t797 - t477 * t796 - t478 * t795 - t479 * t794 + (t504 * t743 + t505 * t742 + t506 * t741 + t507 * t740) * t583, t475 * t993 + t477 * t991 + t478 * t989 + t479 * t987 + (-t504 * t801 - t505 * t800 - t506 * t799 - t507 * t798) * t583, (-t504 * t739 - t505 * t737 - t506 * t735 - t507 * t733) * t583, (-t504 * t738 - t505 * t736 - t506 * t734 - t507 * t732) * t999, (t504 * t615 + t505 * t614 + t506 * t613 + t507 * t612) * t583, (t504 * t659 + t505 * t658 + t506 * t657 + t507 * t656) * t583, (-t504 * t984 - t505 * t982 - t506 * t980 - t507 * t978) * t583, t401 * t865 + t403 * t863 + t404 * t862 + t405 * t861 + (t616 * t680 + t617 * t681 + t618 * t682 + t619 * t683) * t583, t402 * t865 + t406 * t863 + t407 * t862 + t408 * t861 + (-t626 * t680 - t628 * t681 - t629 * t682 - t631 * t683) * t583, 0, -t508, -t509, 0; t398, t415, t394, t395, t414, t411, t412, t413, t436, t452 * t864 + t456 * t860 + t458 * t859 + t460 * t858 + (-t475 * t651 - t477 * t650 - t478 * t649 - t479 * t648) * t583, -t468 * t864 - t470 * t860 - t472 * t859 - t474 * t858 + (t475 * t638 + t477 * t636 + t478 * t634 + t479 * t632) * t583, 0, 0, 0, 0; t476 ^ 2 * t941 + t480 ^ 2 * t927 + t481 ^ 2 * t919 + t482 ^ 2 * t911, (t492 * t788 + t493 * t776 + t494 * t764 + t495 * t752) * t584, (-t747 * t958 - t759 * t959 - t771 * t960 - t783 * t964) * t999, (t748 * t958 + t760 * t959 + t772 * t960 + t784 * t964) * t999, (t492 * t793 + t493 * t781 + t494 * t769 + t495 * t757) * t584, (t492 * t792 + t493 * t780 + t494 * t768 + t495 * t756) * t998, (-t492 * t790 - t493 * t778 - t494 * t766 - t495 * t754) * t998, (-t492 * t785 - t493 * t773 - t494 * t761 - t495 * t749) * t998, (t492 * t516 + t493 * t537 + t494 * t543 + t495 * t549) * t584, -t451 * t864 - t455 * t860 - t457 * t859 - t459 * t858 + (-t476 * t651 - t480 * t650 - t481 * t649 - t482 * t648) * t583, t467 * t864 + t469 * t860 + t471 * t859 + t473 * t858 + (t476 * t638 + t480 * t636 + t481 * t634 + t482 * t632) * t583, 0, 0, 0, t491; t410, t427, t397, t400, t426, t424, t434, t437, 0, ((-t482 * t908 - t503 * t688) * t548 + (-t481 * t916 - t502 * t689) * t542 + (-t480 * t924 - t501 * t690) * t536 + (-t476 * t938 - t500 * t691) * t515) * t583, (t500 * t677 + t501 * t674 + t502 * t671 + t503 * t668) * t583 + t744, 0, 0, 0, 0; t1002, (t430 * t719 + t431 * t718 + t432 * t717 + t433 * t716) * t583, -t476 * t797 - t480 * t796 - t481 * t795 - t482 * t794 + (-t500 * t743 - t501 * t742 - t502 * t741 - t503 * t740) * t583, t476 * t993 + t480 * t991 + t481 * t989 + t482 * t987 + (t500 * t801 + t501 * t800 + t502 * t799 + t503 * t798) * t583, (t500 * t739 + t501 * t737 + t502 * t735 + t503 * t733) * t583, (t500 * t738 + t501 * t736 + t502 * t734 + t503 * t732) * t999, (-t500 * t615 - t501 * t614 - t502 * t613 - t503 * t612) * t583, (-t500 * t659 - t501 * t658 - t502 * t657 - t503 * t656) * t583, (t500 * t984 + t501 * t982 + t502 * t980 + t503 * t978) * t583, t401 * t864 + t403 * t860 + t404 * t859 + t405 * t858 + (-t620 * t680 - t621 * t681 - t622 * t682 - t623 * t683) * t583, t402 * t864 + t406 * t860 + t407 * t859 + t408 * t858 + (t633 * t680 + t635 * t681 + t637 * t682 + t639 * t683) * t583, 0, t509, -t508, 0; t409, t429, t396, t399, t428, t425, t435, t438, 0, -t409 * t583 + t452 * t942 + t456 * t928 + t458 * t920 + t460 * t912, -t468 * t942 - t470 * t928 - t472 * t920 - t474 * t912 + t745, 0, 0, 0, 0; t410, t427, t397, t400, t426, t424, t434, t437, 0, -t410 * t583 - t451 * t942 - t455 * t928 - t457 * t920 - t459 * t912, t467 * t942 + t469 * t928 + t471 * t920 + t473 * t912 + t744, 0, 0, 0, 0; t611, (t911 + t919 + t927 + t941) * t584, (-t515 * t939 - t536 * t925 - t542 * t917 - t548 * t909) * t999, (t827 + t836 + t845 + t856) * t999, (t828 + t837 + t846 + t857) * t584, (t514 * t937 + t525 * t907 + t530 * t901 + t535 * t895) * t998, 0, 0, 0, -0.2e1 * t611 * t583, (t515 * t849 + t536 * t838 + t542 * t829 + t548 * t820) * t999, 0, 0, 0, 1; t1000, (t874 + t875 + t876 + t877) * t583, t1003 * t583 - t430 * t943 - t431 * t929 - t432 * t921 - t433 * t913, t1005 * t583 + t986 + t988 + t990 + t992, (t512 * t877 + t523 * t876 + t528 * t875 + t533 * t874) * t583, (t430 * t513 * t554 + t431 * t524 * t562 + t432 * t529 * t564 + t433 * t534 * t566) * t999, (-t895 * t979 - t901 * t981 - t907 * t983 - t937 * t985) * t583, (-t979 - t981 - t983 - t985) * t583, 0, -t1000 * t583 + t401 * t942 + t403 * t928 + t404 * t920 + t405 * t912, t1004 * t583 + t402 * t942 + t406 * t928 + t407 * t920 + t408 * t912, 0, 0, 0, 0; t1001, (-t443 * t697 - t444 * t712 - t445 * t707 - t446 * t702) * t584, ((-t479 * t872 + t507 * t878) * t890 + (-t478 * t866 + t506 * t879) * t896 + (-t477 * t868 + t505 * t880) * t902 + (-t475 * t870 + t504 * t881) * t932) * t583, ((t479 * t873 - t680 * t947) * t883 + (t478 * t867 - t681 * t949) * t885 + (t477 * t869 - t682 * t951) * t887 + (t475 * t871 - t683 * t953) * t889) * t583, (-t443 * t700 - t444 * t715 - t445 * t710 - t446 * t705) * t584, (-t443 * t699 - t444 * t714 - t445 * t709 - t446 * t704) * t998, (t443 * t698 + t444 * t713 + t445 * t708 + t446 * t703) * t998, (t443 * t696 + t444 * t711 + t445 * t706 + t446 * t701) * t998, (-t443 * t947 - t444 * t953 - t445 * t951 - t446 * t949) * t584, t683 * t452 + t682 * t456 + t681 * t458 + t680 * t460 + (-t475 * t654 - t477 * t653 - t478 * t652 - t479 * t655) * t583, -t683 * t468 - t682 * t470 - t681 * t472 - t680 * t474 + (t475 * t642 + t477 * t641 + t478 * t640 + t479 * t643) * t583, 0, -t508, -t509, 0; t1002, (t512 * t728 + t523 * t725 + t528 * t722 + t533 * t731) * t584, ((-t482 * t872 - t503 * t878) * t890 + (-t481 * t866 - t502 * t879) * t896 + (-t480 * t868 - t501 * t880) * t902 + (-t476 * t870 - t500 * t881) * t932) * t583, ((t482 * t873 + t680 * t954) * t883 + (t481 * t867 + t681 * t955) * t885 + (t480 * t869 + t682 * t956) * t887 + (t476 * t871 + t683 * t957) * t889) * t583, (t510 * t728 + t521 * t725 + t526 * t722 + t531 * t731) * t584, (t756 * t977 + t768 * t968 + t780 * t971 + t792 * t974) * t998, (-t754 * t977 - t766 * t968 - t778 * t971 - t790 * t974) * t998, (-t443 * t716 - t444 * t719 - t445 * t718 - t446 * t717) * t998, (t443 * t954 + t444 * t957 + t445 * t956 + t446 * t955) * t584, -t683 * t451 - t682 * t455 - t681 * t457 - t680 * t459 + (-t476 * t654 - t480 * t653 - t481 * t652 - t482 * t655) * t583, t683 * t467 + t682 * t469 + t681 * t471 + t680 * t473 + (t476 * t642 + t480 * t641 + t481 * t640 + t482 * t643) * t583, 0, t509, -t508, 0; t1000, (t443 * t746 + t444 * t782 + t445 * t770 + t446 * t758) * t584, (-t443 * t566 * t823 - t444 * t554 * t852 - t445 * t562 * t841 - t446 * t564 * t832 + t1003) * t583, (t1005 + t721 + t724 + t727 + t730) * t583, (t443 * t755 + t444 * t791 + t445 * t779 + t446 * t767) * t584, (t443 * t753 + t444 * t789 + t445 * t777 + t446 * t765) * t998, (-t443 * t826 - t444 * t855 - t445 * t844 - t446 * t835) * t584, (-t443 * t915 - t444 * t945 - t445 * t931 - t446 * t923) * t584, 0, (-t688 * t976 - t689 * t967 - t690 * t970 - t691 * t973 - t1000) * t583, (t443 * t668 + t444 * t677 + t445 * t674 + t446 * t671 + t1004) * t583, 0, 0, 0, 0; t680 ^ 2 + t681 ^ 2 + t682 ^ 2 + t683 ^ 2, (t430 * t726 + t431 * t723 + t432 * t720 + t433 * t729) * t583, -t683 * t992 - t682 * t990 - t681 * t988 - t680 * t986 + (-t680 * t730 - t681 * t721 - t682 * t724 - t683 * t727) * t583, t683 * t430 * t555 + t682 * t431 * t563 + t681 * t432 * t565 + t680 * t433 * t567 + (t443 * t798 + t444 * t801 + t445 * t800 + t446 * t799) * t583, (t443 * t733 + t444 * t739 + t445 * t737 + t446 * t735) * t583, (t443 * t732 + t444 * t738 + t445 * t736 + t446 * t734) * t999, (-t443 * t612 - t444 * t615 - t445 * t614 - t446 * t613) * t583, (-t443 * t656 - t444 * t659 - t445 * t658 - t446 * t657) * t583, (t439 * t973 + t440 * t970 + t441 * t967 + t442 * t976) * t583, t683 * t401 + t682 * t403 + t681 * t404 + t680 * t405 + (-t660 * t994 - t661 * t995 - t662 * t996 - t663 * t997) * t583, t683 * t402 + t682 * t406 + t681 * t407 + t680 * t408 + (t670 * t994 + t673 * t995 + t676 * t996 + t679 * t997) * t583, 1, 0, 0, 0;];
tau_reg  = t1;