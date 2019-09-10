% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x10]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRR1G1P1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:34
% EndTime: 2019-05-03 15:38:37
% DurationCPUTime: 3.70s
% Computational Cost: add. (24372->237), mult. (38819->488), div. (2856->9), fcn. (29350->26), ass. (0->226)
t663 = xP(3);
t639 = sin(t663);
t640 = cos(t663);
t664 = koppelP(3,2);
t667 = koppelP(3,1);
t609 = t639 * t667 + t640 * t664;
t660 = xDP(3);
t662 = xDP(1);
t582 = t609 * t660 - t662;
t612 = -t639 * t664 + t640 * t667;
t661 = xDP(2);
t585 = t612 * t660 + t661;
t642 = qJ(1,3) + qJ(2,3);
t627 = sin(t642);
t630 = cos(t642);
t645 = legFrame(3,3);
t633 = sin(t645);
t636 = cos(t645);
t791 = (t582 * t636 - t585 * t633) * t630 - (t582 * t633 + t585 * t636) * t627;
t665 = koppelP(2,2);
t668 = koppelP(2,1);
t610 = t639 * t668 + t640 * t665;
t583 = t610 * t660 - t662;
t613 = -t639 * t665 + t640 * t668;
t586 = t613 * t660 + t661;
t643 = qJ(1,2) + qJ(2,2);
t628 = sin(t643);
t631 = cos(t643);
t646 = legFrame(2,3);
t634 = sin(t646);
t637 = cos(t646);
t790 = (t583 * t637 - t586 * t634) * t631 - (t583 * t634 + t586 * t637) * t628;
t666 = koppelP(1,2);
t669 = koppelP(1,1);
t611 = t639 * t669 + t640 * t666;
t584 = t611 * t660 - t662;
t614 = -t639 * t666 + t640 * t669;
t587 = t614 * t660 + t661;
t644 = qJ(1,1) + qJ(2,1);
t629 = sin(t644);
t632 = cos(t644);
t647 = legFrame(1,3);
t635 = sin(t647);
t638 = cos(t647);
t789 = (t584 * t638 - t587 * t635) * t632 - (t584 * t635 + t587 * t638) * t629;
t649 = sin(qJ(1,3));
t655 = cos(qJ(1,3));
t615 = -t667 * t649 + t655 * t664;
t618 = t649 * t664 + t655 * t667;
t788 = (t615 * t640 + t618 * t639) * t636 - t633 * (-t615 * t639 + t618 * t640);
t651 = sin(qJ(1,2));
t657 = cos(qJ(1,2));
t616 = -t668 * t651 + t657 * t665;
t619 = t651 * t665 + t657 * t668;
t787 = (t616 * t640 + t619 * t639) * t637 - t634 * (-t616 * t639 + t619 * t640);
t653 = sin(qJ(1,1));
t659 = cos(qJ(1,1));
t617 = -t669 * t653 + t659 * t666;
t620 = t653 * t666 + t659 * t669;
t786 = (t617 * t640 + t620 * t639) * t638 - t635 * (-t617 * t639 + t620 * t640);
t606 = t627 * t655 - t630 * t649;
t785 = 0.1e1 / t606;
t607 = t628 * t657 - t631 * t651;
t784 = 0.1e1 / t607;
t608 = t629 * t659 - t632 * t653;
t783 = 0.1e1 / t608;
t671 = 0.1e1 / pkin(2);
t549 = pkin(1) * ((-t649 * t661 - t655 * t662) * t636 - t633 * (-t649 * t662 + t655 * t661) + t788 * t660) + t791 * pkin(2);
t672 = 0.1e1 / pkin(1);
t752 = t785 * t672;
t727 = t549 * t752;
t546 = t671 * t727;
t555 = t791 * t752;
t541 = -t555 + t546;
t691 = t627 * t649 + t630 * t655;
t770 = t791 * t785;
t534 = -pkin(2) * t541 + t691 * t770;
t598 = 0.1e1 / t606 ^ 2;
t673 = 0.1e1 / pkin(1) ^ 2;
t588 = t627 * t636 + t630 * t633;
t570 = pkin(1) * (t633 * t655 + t636 * t649) + t588 * pkin(2);
t589 = -t627 * t633 + t630 * t636;
t573 = pkin(1) * (-t633 * t649 + t636 * t655) + t589 * pkin(2);
t676 = (-t570 * t609 - t573 * t612) * t671 * t785;
t697 = -t588 * t609 - t589 * t612;
t685 = t697 * t785;
t773 = t541 * t549;
t730 = t785 * t773;
t718 = t785 * t730;
t719 = (pkin(1) * t691 + pkin(2)) * t671 * t773;
t641 = t660 ^ 2;
t744 = t641 * t672;
t670 = pkin(2) ^ 2;
t743 = 0.2e1 * pkin(2);
t776 = (t541 * t670 + (-t770 + t691 * (-t555 + t546 / 0.2e1) * t743) * pkin(1)) * t671;
t522 = (0.2e1 * t718 + (-t719 - (-0.2e1 * t534 - t776) * t791) * t598) * t673 + (-t676 + 0.2e1 * t685) * t744;
t782 = t522 * t785;
t550 = pkin(1) * ((-t651 * t661 - t657 * t662) * t637 - t634 * (-t651 * t662 + t657 * t661) + t787 * t660) + t790 * pkin(2);
t750 = t784 * t672;
t726 = t550 * t750;
t547 = t671 * t726;
t556 = t790 * t750;
t543 = -t556 + t547;
t690 = t628 * t651 + t631 * t657;
t769 = t790 * t784;
t535 = -pkin(2) * t543 + t690 * t769;
t600 = 0.1e1 / t607 ^ 2;
t590 = t628 * t637 + t631 * t634;
t571 = pkin(1) * (t634 * t657 + t637 * t651) + t590 * pkin(2);
t591 = -t628 * t634 + t631 * t637;
t574 = pkin(1) * (-t634 * t651 + t637 * t657) + t591 * pkin(2);
t675 = (-t571 * t610 - t574 * t613) * t671 * t784;
t696 = -t590 * t610 - t591 * t613;
t684 = t696 * t784;
t772 = t543 * t550;
t729 = t784 * t772;
t714 = t784 * t729;
t715 = (pkin(1) * t690 + pkin(2)) * t671 * t772;
t775 = (t543 * t670 + (-t769 + t690 * (-t556 + t547 / 0.2e1) * t743) * pkin(1)) * t671;
t523 = (0.2e1 * t714 + (-t715 - (-0.2e1 * t535 - t775) * t790) * t600) * t673 + (-t675 + 0.2e1 * t684) * t744;
t781 = t523 * t784;
t551 = pkin(1) * ((-t653 * t661 - t659 * t662) * t638 - t635 * (-t653 * t662 + t659 * t661) + t786 * t660) + t789 * pkin(2);
t748 = t783 * t672;
t725 = t551 * t748;
t548 = t671 * t725;
t557 = t789 * t748;
t545 = -t557 + t548;
t689 = t629 * t653 + t632 * t659;
t768 = t789 * t783;
t536 = -pkin(2) * t545 + t689 * t768;
t602 = 0.1e1 / t608 ^ 2;
t592 = t629 * t638 + t632 * t635;
t572 = pkin(1) * (t635 * t659 + t638 * t653) + t592 * pkin(2);
t593 = -t629 * t635 + t632 * t638;
t575 = pkin(1) * (-t635 * t653 + t638 * t659) + t593 * pkin(2);
t674 = (-t572 * t611 - t575 * t614) * t671 * t783;
t695 = -t592 * t611 - t593 * t614;
t683 = t695 * t783;
t771 = t545 * t551;
t728 = t783 * t771;
t710 = t783 * t728;
t711 = (pkin(1) * t689 + pkin(2)) * t671 * t771;
t774 = (t545 * t670 + (-t768 + t689 * (-t557 + t548 / 0.2e1) * t743) * pkin(1)) * t671;
t524 = (0.2e1 * t710 + (-t711 - (-0.2e1 * t536 - t774) * t789) * t602) * t673 + (-t674 + 0.2e1 * t683) * t744;
t780 = t524 * t783;
t525 = (t718 + (-t719 - (-t534 - t776) * t791) * t598) * t673 + (-t676 + t685) * t744;
t779 = t525 * t785;
t526 = (t714 + (-t715 - (-t535 - t775) * t790) * t600) * t673 + (-t675 + t684) * t744;
t778 = t526 * t784;
t527 = (t710 + (-t711 - (-t536 - t774) * t789) * t602) * t673 + (-t674 + t683) * t744;
t777 = t527 * t783;
t564 = (t609 * t636 - t612 * t633) * t630 - (t609 * t633 + t612 * t636) * t627;
t767 = t564 * t785;
t565 = (t610 * t637 - t613 * t634) * t631 - (t610 * t634 + t613 * t637) * t628;
t766 = t565 * t784;
t566 = (t611 * t638 - t614 * t635) * t632 - (t611 * t635 + t614 * t638) * t629;
t765 = t566 * t783;
t761 = t570 * t785;
t760 = t571 * t784;
t759 = t572 * t783;
t758 = t588 * t785;
t757 = t589 * t785;
t756 = t590 * t784;
t755 = t591 * t784;
t754 = t592 * t783;
t753 = t593 * t783;
t751 = t598 * t673;
t749 = t600 * t673;
t747 = t602 * t673;
t746 = t639 * t641;
t745 = t640 * t641;
t742 = t522 * t758;
t648 = sin(qJ(2,3));
t741 = t648 * t782;
t654 = cos(qJ(2,3));
t740 = t654 * t782;
t739 = t523 * t756;
t650 = sin(qJ(2,2));
t738 = t650 * t781;
t656 = cos(qJ(2,2));
t737 = t656 * t781;
t736 = t524 * t754;
t652 = sin(qJ(2,1));
t735 = t652 * t780;
t658 = cos(qJ(2,1));
t734 = t658 * t780;
t528 = t534 * t791 * t751 + (t673 * t730 + t697 * t744) * t785;
t733 = t528 * t761;
t529 = t535 * t790 * t749 + (t673 * t729 + t696 * t744) * t784;
t732 = t529 * t760;
t530 = t536 * t789 * t747 + (t673 * t728 + t695 * t744) * t783;
t731 = t530 * t759;
t558 = t791 ^ 2;
t724 = t558 * t751;
t559 = t790 ^ 2;
t723 = t559 * t749;
t560 = t789 ^ 2;
t722 = t560 * t747;
t540 = -0.2e1 * t555 + t546;
t721 = t540 * t549 * t588 * t598;
t720 = t540 * t727;
t542 = -0.2e1 * t556 + t547;
t717 = t542 * t550 * t590 * t600;
t716 = t542 * t726;
t544 = -0.2e1 * t557 + t548;
t713 = t544 * t551 * t592 * t602;
t712 = t544 * t725;
t709 = t558 * t598 * t761;
t708 = t559 * t600 * t760;
t707 = t560 * t602 * t759;
t706 = t648 * t720;
t705 = t654 * t720;
t704 = t650 * t716;
t703 = t656 * t716;
t702 = t652 * t712;
t701 = t658 * t712;
t682 = -t528 * t648 + t654 * t724;
t681 = t528 * t654 + t648 * t724;
t680 = -t529 * t650 + t656 * t723;
t679 = t529 * t656 + t650 * t723;
t678 = -t530 * t652 + t658 * t722;
t677 = t530 * t658 + t652 * t722;
t554 = -pkin(1) * t786 - t566 * pkin(2);
t553 = -pkin(1) * t787 - t565 * pkin(2);
t552 = -pkin(1) * t788 - t564 * pkin(2);
t1 = [(t528 * t757 + t529 * t755 + t530 * t753) * t672, 0, 0, (t525 * t757 + t526 * t755 + t527 * t753 + (-t573 * t779 - t574 * t778 - t575 * t777) * t671) * t672, t589 * t740 + t591 * t737 + t593 * t734 + (-(t575 * t677 + t593 * t702) * t783 - (t574 * t679 + t591 * t704) * t784 - (t573 * t681 + t589 * t706) * t785) * t671, -t589 * t741 - t591 * t738 - t593 * t735 + (-(t575 * t678 + t593 * t701) * t783 - (t574 * t680 + t591 * t703) * t784 - (t573 * t682 + t589 * t705) * t785) * t671, 0, -t745, t746, 0; (t528 * t758 + t529 * t756 + t530 * t754) * t672, 0, 0, (t525 * t758 + t526 * t756 + t527 * t754 + (-t525 * t761 - t526 * t760 - t527 * t759) * t671) * t672, t654 * t742 + t656 * t739 + t658 * t736 + (-t654 * t733 - t656 * t732 - t658 * t731 + (-t648 * t709 - t650 * t708 - t652 * t707) * t673 + (-t648 * t721 - t650 * t717 - t652 * t713) * t672) * t671, -t648 * t742 - t650 * t739 - t652 * t736 + (t648 * t733 + t650 * t732 + t652 * t731 + (-t654 * t709 - t656 * t708 - t658 * t707) * t673 + (-t654 * t721 - t656 * t717 - t658 * t713) * t672) * t671, 0, -t746, -t745, 0; (-t528 * t767 - t529 * t766 - t530 * t765) * t672, 0, 0, (-t525 * t767 - t526 * t766 - t527 * t765 + (-t552 * t779 - t553 * t778 - t554 * t777) * t671) * t672, -t564 * t740 - t565 * t737 - t566 * t734 + (-(t554 * t677 - t566 * t702) * t783 - (t553 * t679 - t565 * t704) * t784 - (t552 * t681 - t564 * t706) * t785) * t671, t564 * t741 + t565 * t738 + t566 * t735 + (-(t554 * t678 - t566 * t701) * t783 - (t553 * t680 - t565 * t703) * t784 - (t552 * t682 - t564 * t705) * t785) * t671, 0, 0, 0, 0;];
tau_reg  = t1;
