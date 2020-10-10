% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRRR1G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:12
% EndTime: 2020-03-09 20:34:13
% DurationCPUTime: 1.83s
% Computational Cost: add. (1899->210), mult. (5679->501), div. (1878->17), fcn. (5343->18), ass. (0->222)
t617 = legFrame(3,3);
t584 = sin(t617);
t587 = cos(t617);
t632 = xDP(2);
t633 = xDP(1);
t578 = t584 * t633 - t587 * t632;
t626 = cos(qJ(3,3));
t620 = sin(qJ(3,3));
t627 = cos(qJ(2,3));
t734 = t620 * t627;
t557 = (-t584 * t632 - t587 * t633) * t626 - t578 * t734;
t621 = sin(qJ(2,3));
t591 = 0.1e1 / t621;
t592 = 0.1e1 / t621 ^ 2;
t602 = 0.1e1 / t626;
t645 = t626 ^ 2;
t603 = 0.1e1 / t645;
t604 = t602 * t603;
t749 = t591 * t627;
t707 = t557 * t749;
t735 = t620 * t621;
t634 = 0.1e1 / pkin(2);
t779 = t634 ^ 2;
t521 = (-(-t578 * t735 + t707) * t592 * t557 - (-t557 * t620 + t578 * t627) * t591 * t578) * t604 * t602 * t779;
t605 = 0.1e1 / t645 ^ 2;
t554 = t557 ^ 2;
t635 = 0.1e1 / pkin(2) ^ 2;
t712 = t554 * t592 * t635;
t548 = t605 * t712;
t575 = t578 ^ 2;
t755 = t575 * t603;
t704 = t635 * t755;
t542 = t548 + t704;
t778 = -0.2e1 * t635;
t662 = t578 * t707 * t778;
t729 = t626 * t627;
t782 = -t542 * t621 * t626 + t620 * t604 * t662 + t521 * t729;
t618 = legFrame(2,3);
t585 = sin(t618);
t588 = cos(t618);
t579 = t585 * t633 - t588 * t632;
t628 = cos(qJ(3,2));
t622 = sin(qJ(3,2));
t629 = cos(qJ(2,2));
t732 = t622 * t629;
t558 = (-t585 * t632 - t588 * t633) * t628 - t579 * t732;
t623 = sin(qJ(2,2));
t595 = 0.1e1 / t623;
t596 = 0.1e1 / t623 ^ 2;
t607 = 0.1e1 / t628;
t649 = t628 ^ 2;
t608 = 0.1e1 / t649;
t609 = t607 * t608;
t746 = t595 * t629;
t706 = t558 * t746;
t733 = t622 * t623;
t522 = (-(-t579 * t733 + t706) * t596 * t558 - (-t558 * t622 + t579 * t629) * t595 * t579) * t609 * t607 * t779;
t610 = 0.1e1 / t649 ^ 2;
t555 = t558 ^ 2;
t710 = t555 * t596 * t635;
t549 = t610 * t710;
t576 = t579 ^ 2;
t753 = t576 * t608;
t702 = t635 * t753;
t543 = t549 + t702;
t661 = t579 * t706 * t778;
t728 = t628 * t629;
t781 = -t543 * t623 * t628 + t622 * t609 * t661 + t522 * t728;
t619 = legFrame(1,3);
t586 = sin(t619);
t589 = cos(t619);
t580 = t586 * t633 - t589 * t632;
t630 = cos(qJ(3,1));
t624 = sin(qJ(3,1));
t631 = cos(qJ(2,1));
t730 = t624 * t631;
t559 = (-t586 * t632 - t589 * t633) * t630 - t580 * t730;
t625 = sin(qJ(2,1));
t599 = 0.1e1 / t625;
t600 = 0.1e1 / t625 ^ 2;
t612 = 0.1e1 / t630;
t653 = t630 ^ 2;
t613 = 0.1e1 / t653;
t614 = t612 * t613;
t743 = t599 * t631;
t705 = t559 * t743;
t731 = t624 * t625;
t523 = (-(-t580 * t731 + t705) * t600 * t559 - (-t559 * t624 + t580 * t631) * t599 * t580) * t614 * t612 * t779;
t615 = 0.1e1 / t653 ^ 2;
t556 = t559 ^ 2;
t708 = t556 * t600 * t635;
t550 = t615 * t708;
t577 = t580 ^ 2;
t751 = t577 * t613;
t700 = t635 * t751;
t544 = t550 + t700;
t660 = t580 * t705 * t778;
t727 = t630 * t631;
t780 = -t544 * t625 * t630 + t624 * t614 * t660 + t523 * t727;
t777 = t521 * t602;
t776 = t521 * t620;
t775 = t521 * t627;
t774 = t522 * t607;
t773 = t522 * t622;
t772 = t522 * t629;
t771 = t523 * t612;
t770 = t523 * t624;
t769 = t523 * t631;
t593 = t591 * t592;
t756 = t575 * t591;
t536 = (t554 * t593 + t756) * t604 * t634;
t768 = t536 * t602;
t767 = t536 * t603;
t597 = t595 * t596;
t754 = t576 * t595;
t537 = (t555 * t597 + t754) * t609 * t634;
t766 = t537 * t607;
t765 = t537 * t608;
t601 = t599 * t600;
t752 = t577 * t599;
t538 = (t556 * t601 + t752) * t614 * t634;
t764 = t538 * t612;
t763 = t538 * t613;
t762 = t554 * t627;
t761 = t555 * t629;
t760 = t556 * t631;
t759 = t557 * t578;
t758 = t558 * t579;
t757 = t559 * t580;
t750 = t591 * t602;
t606 = t602 * t605;
t748 = t592 * t606;
t747 = t595 * t607;
t611 = t607 * t610;
t745 = t596 * t611;
t744 = t599 * t612;
t616 = t612 * t615;
t742 = t600 * t616;
t740 = t605 * t620;
t738 = t610 * t622;
t736 = t615 * t624;
t527 = t542 * t735 + t603 * t662;
t528 = t543 * t733 + t608 * t661;
t529 = t544 * t731 + t613 * t660;
t636 = t634 * t635;
t726 = 0.2e1 * t636;
t725 = 0.2e1 * t759;
t724 = 0.2e1 * t758;
t723 = 0.2e1 * t757;
t722 = t521 * t591 * t603;
t721 = t602 * t776;
t720 = t522 * t595 * t608;
t719 = t607 * t773;
t718 = t523 * t599 * t613;
t717 = t612 * t770;
t716 = t536 * t750;
t715 = t537 * t747;
t714 = t538 * t744;
t713 = t554 * t748;
t711 = t555 * t745;
t709 = t556 * t742;
t703 = t575 * t740;
t701 = t576 * t738;
t699 = t577 * t736;
t698 = t603 * t749;
t697 = t592 * t740;
t696 = t608 * t746;
t695 = t596 * t738;
t694 = t613 * t743;
t693 = t600 * t736;
t692 = 0.2e1 * t591 * t776;
t691 = 0.2e1 * t595 * t773;
t690 = 0.2e1 * t599 * t770;
t590 = t620 ^ 2;
t689 = t590 * t722;
t688 = t749 * t777;
t594 = t622 ^ 2;
t687 = t594 * t720;
t686 = t746 * t774;
t598 = t624 ^ 2;
t685 = t598 * t718;
t684 = t743 * t771;
t683 = t536 * t698;
t682 = t537 * t696;
t681 = t538 * t694;
t680 = t593 * t606 * t762;
t679 = t597 * t611 * t761;
t678 = t601 * t616 * t760;
t677 = t575 * t590 * t604 * t621;
t676 = t576 * t594 * t609 * t623;
t675 = t577 * t598 * t614 * t625;
t674 = t620 * t698;
t673 = t622 * t696;
t672 = t624 * t694;
t671 = (-t635 * t676 + t781) * t747;
t670 = (-t635 * t675 + t780) * t744;
t669 = (-t635 * t677 + t782) * t750;
t668 = ((-t621 * t704 - t775) * t620 + t527) * t750;
t667 = ((-t623 * t702 - t772) * t622 + t528) * t747;
t666 = ((-t625 * t700 - t769) * t624 + t529) * t744;
t665 = (-0.1e1 + 0.2e1 * t645) * t748 * t759;
t664 = (-0.1e1 + 0.2e1 * t649) * t745 * t758;
t663 = (-0.1e1 + 0.2e1 * t653) * t742 * t757;
t659 = (t590 * t606 + t604) * t756;
t658 = (t594 * t611 + t609) * t754;
t657 = (t598 * t616 + t614) * t752;
t574 = t586 * t727 - t589 * t624;
t573 = t585 * t728 - t588 * t622;
t572 = t584 * t729 - t587 * t620;
t571 = -t586 * t630 + t589 * t730;
t570 = t586 * t624 + t589 * t727;
t569 = -t585 * t628 + t588 * t732;
t568 = t585 * t622 + t588 * t728;
t567 = -t584 * t626 + t587 * t734;
t566 = t584 * t620 + t587 * t729;
t565 = -t586 * t730 - t589 * t630;
t564 = -t585 * t732 - t588 * t628;
t563 = -t584 * t734 - t587 * t626;
t532 = -0.2e1 * t613 * t708 + t550;
t531 = -0.2e1 * t608 * t710 + t549;
t530 = -0.2e1 * t603 * t712 + t548;
t1 = [t566 * t716 + t568 * t715 + t570 * t714, (t563 * t722 + t564 * t720 + t565 * t718) * t634, t566 * t688 + t568 * t686 + t570 * t684 + (-t566 * t713 - t568 * t711 - t570 * t709) * t635 + (t563 * t683 + t564 * t682 + t565 * t681) * t634, -t566 * t777 - t568 * t774 - t570 * t771 + (-t566 * t680 - t568 * t679 - t570 * t678) * t635 + (-t563 * t767 - t564 * t765 - t565 * t763) * t634, (t563 * t689 + t564 * t687 + t565 * t685) * t634 + ((-t556 * t586 + t565 * t723) * t693 + (-t555 * t585 + t564 * t724) * t695 + (-t554 * t584 + t563 * t725) * t697) * t636, (t563 * t665 + t564 * t664 + t565 * t663) * t726 + ((t532 * t586 + t565 * t690) * t612 + (t531 * t585 + t564 * t691) * t607 + (t530 * t584 + t563 * t692) * t602) * t634, (t584 * t721 + t585 * t719 + t586 * t717) * t634 + (t563 * t659 + t564 * t658 + t565 * t657) * t636, (t521 * t584 + t522 * t585 + t523 * t586) * t634, (t584 * t703 + t585 * t701 + t586 * t699) * t636, t570 * t670 + t568 * t671 + t566 * t669 + ((t565 * t743 - t586 * t731) * t764 + (t564 * t746 - t585 * t733) * t766 + (t563 * t749 - t584 * t735) * t768) * t634, t570 * t666 + t568 * t667 + t566 * t668 + ((-t565 * t672 - t586 * t625) * t538 + (-t564 * t673 - t585 * t623) * t537 + (-t563 * t674 - t584 * t621) * t536) * t634, 0; t572 * t716 + t573 * t715 + t574 * t714, (t567 * t722 + t569 * t720 + t571 * t718) * t634, t572 * t688 + t573 * t686 + t574 * t684 + (-t572 * t713 - t573 * t711 - t574 * t709) * t635 + (t567 * t683 + t569 * t682 + t571 * t681) * t634, -t572 * t777 - t573 * t774 - t574 * t771 + (-t572 * t680 - t573 * t679 - t574 * t678) * t635 + (-t567 * t767 - t569 * t765 - t571 * t763) * t634, (t567 * t689 + t569 * t687 + t571 * t685) * t634 + ((t556 * t589 + t571 * t723) * t693 + (t555 * t588 + t569 * t724) * t695 + (t554 * t587 + t567 * t725) * t697) * t636, (t567 * t665 + t569 * t664 + t571 * t663) * t726 + ((-t532 * t589 + t571 * t690) * t612 + (-t531 * t588 + t569 * t691) * t607 + (-t530 * t587 + t567 * t692) * t602) * t634, (-t587 * t721 - t588 * t719 - t589 * t717) * t634 + (t567 * t659 + t569 * t658 + t571 * t657) * t636, (-t521 * t587 - t522 * t588 - t523 * t589) * t634, (-t587 * t703 - t588 * t701 - t589 * t699) * t636, t574 * t670 + t573 * t671 + t572 * t669 + ((t571 * t743 + t589 * t731) * t764 + (t569 * t746 + t588 * t733) * t766 + (t567 * t749 + t587 * t735) * t768) * t634, t574 * t666 + t573 * t667 + t572 * t668 + ((-t571 * t672 + t589 * t625) * t538 + (-t569 * t673 + t588 * t623) * t537 + (-t567 * t674 + t587 * t621) * t536) * t634, 0; t538 + t537 + t536, 0, t775 + t772 + t769 + (-t554 * t591 * t605 - t555 * t595 * t610 - t556 * t599 * t615) * t635, -t521 * t621 - t522 * t623 - t523 * t625 + (-t592 * t605 * t762 - t596 * t610 * t761 - t600 * t615 * t760) * t635, 0, 0, 0, 0, 0, (-t675 - t676 - t677) * t635 + t780 + t781 + t782, -t521 * t734 - t522 * t732 - t523 * t730 + (-t731 * t751 - t733 * t753 - t735 * t755) * t635 + t529 + t528 + t527, 0;];
tau_reg  = t1;
