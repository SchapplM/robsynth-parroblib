% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR1G3P3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G3P3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:57
% EndTime: 2020-03-09 21:06:58
% DurationCPUTime: 1.73s
% Computational Cost: add. (19011->153), mult. (11601->332), div. (2418->9), fcn. (11529->24), ass. (0->173)
t599 = pkin(7) + qJ(2,3);
t590 = qJ(3,3) + t599;
t578 = sin(t590);
t581 = cos(t590);
t611 = xDP(3);
t602 = legFrame(3,2);
t593 = sin(t602);
t596 = cos(t602);
t612 = xDP(2);
t613 = xDP(1);
t632 = t593 * t612 - t596 * t613;
t542 = t578 * t611 + t632 * t581;
t600 = pkin(7) + qJ(2,2);
t591 = qJ(3,2) + t600;
t579 = sin(t591);
t582 = cos(t591);
t603 = legFrame(2,2);
t594 = sin(t603);
t597 = cos(t603);
t631 = t594 * t612 - t597 * t613;
t543 = t579 * t611 + t631 * t582;
t601 = pkin(7) + qJ(2,1);
t592 = qJ(3,1) + t601;
t580 = sin(t592);
t583 = cos(t592);
t604 = legFrame(1,2);
t595 = sin(t604);
t598 = cos(t604);
t630 = t595 * t612 - t598 * t613;
t544 = t580 * t611 + t630 * t583;
t616 = 0.1e1 / pkin(2);
t722 = t616 ^ 2;
t584 = sin(t599);
t587 = cos(t599);
t560 = t578 * t587 - t584 * t581;
t718 = 0.1e1 / t560;
t721 = t596 * t718;
t585 = sin(t600);
t588 = cos(t600);
t561 = t579 * t588 - t585 * t582;
t717 = 0.1e1 / t561;
t720 = t597 * t717;
t586 = sin(t601);
t589 = cos(t601);
t562 = t580 * t589 - t586 * t583;
t716 = 0.1e1 / t562;
t719 = t598 * t716;
t686 = t716 * t595;
t687 = t717 * t594;
t688 = t718 * t593;
t715 = 0.1e1 / t560 ^ 2;
t714 = 0.1e1 / t561 ^ 2;
t713 = 0.1e1 / t562 ^ 2;
t615 = 0.1e1 / pkin(3);
t536 = -(-t584 * t611 - t587 * t632) * pkin(2) + t542 * pkin(3);
t700 = t536 * t718;
t675 = t615 * t700;
t694 = t542 * t718;
t531 = (t675 - t694) * t616;
t635 = t578 * t584 + t581 * t587;
t521 = t531 * pkin(3) - t635 * t694;
t655 = t531 * t718 * t700;
t617 = 0.1e1 / pkin(2) ^ 2;
t679 = t615 * t617;
t524 = (t635 * pkin(2) + pkin(3)) * t655 * t679;
t693 = t542 * t715;
t614 = pkin(3) ^ 2;
t657 = 0.2e1 * pkin(3) * t616;
t703 = (-t531 * t614 + (t694 - t635 * (-t694 + t675 / 0.2e1) * t657) * pkin(2)) * t615;
t509 = -t524 + (t655 + (-t521 - t703) * t693) * t617;
t572 = pkin(2) * t587 + pkin(3) * t581;
t712 = t509 * t572;
t537 = -(-t585 * t611 - t588 * t631) * pkin(2) + t543 * pkin(3);
t699 = t537 * t717;
t674 = t615 * t699;
t692 = t543 * t717;
t533 = (t674 - t692) * t616;
t634 = t579 * t585 + t582 * t588;
t522 = t533 * pkin(3) - t634 * t692;
t653 = t533 * t717 * t699;
t525 = (t634 * pkin(2) + pkin(3)) * t653 * t679;
t691 = t543 * t714;
t702 = (-t533 * t614 + (t692 - t634 * (-t692 + t674 / 0.2e1) * t657) * pkin(2)) * t615;
t511 = -t525 + (t653 + (-t522 - t702) * t691) * t617;
t573 = pkin(2) * t588 + pkin(3) * t582;
t711 = t511 * t573;
t538 = -(-t586 * t611 - t589 * t630) * pkin(2) + t544 * pkin(3);
t698 = t538 * t716;
t673 = t615 * t698;
t690 = t544 * t716;
t535 = (t673 - t690) * t616;
t633 = t580 * t586 + t583 * t589;
t523 = t535 * pkin(3) - t633 * t690;
t651 = t535 * t716 * t698;
t526 = (t633 * pkin(2) + pkin(3)) * t651 * t679;
t689 = t544 * t713;
t701 = (-t535 * t614 + (t690 - t633 * (-t690 + t673 / 0.2e1) * t657) * pkin(2)) * t615;
t513 = -t526 + (t651 + (-t523 - t701) * t689) * t617;
t574 = pkin(2) * t589 + pkin(3) * t583;
t710 = t513 * t574;
t515 = (-t521 * t693 + t655) * t617;
t605 = sin(qJ(3,3));
t709 = t515 * t605;
t608 = cos(qJ(3,3));
t708 = t515 * t608;
t516 = (-t522 * t691 + t653) * t617;
t606 = sin(qJ(3,2));
t707 = t516 * t606;
t609 = cos(qJ(3,2));
t706 = t516 * t609;
t517 = (-t523 * t689 + t651) * t617;
t607 = sin(qJ(3,1));
t705 = t517 * t607;
t610 = cos(qJ(3,1));
t704 = t517 * t610;
t697 = t542 ^ 2 * t617;
t696 = t543 ^ 2 * t617;
t695 = t544 ^ 2 * t617;
t685 = t718 * t578;
t684 = t717 * t579;
t683 = t716 * t580;
t678 = (t675 - 0.2e1 * t694) * t722 * t536;
t677 = (t674 - 0.2e1 * t692) * t722 * t537;
t676 = (t673 - 0.2e1 * t690) * t722 * t538;
t672 = t605 * t697;
t671 = t608 * t697;
t670 = t606 * t696;
t669 = t609 * t696;
t668 = t607 * t695;
t667 = t610 * t695;
t666 = t581 * t721;
t665 = t582 * t720;
t664 = t583 * t719;
t663 = t581 * t688;
t662 = t582 * t687;
t661 = t583 * t686;
t510 = -t524 + (0.2e1 * t655 + (-0.2e1 * t521 - t703) * t693) * t617;
t660 = t510 * t685;
t512 = -t525 + (0.2e1 * t653 + (-0.2e1 * t522 - t702) * t691) * t617;
t659 = t512 * t684;
t514 = -t526 + (0.2e1 * t651 + (-0.2e1 * t523 - t701) * t689) * t617;
t658 = t514 * t683;
t656 = t581 * t678;
t654 = t582 * t677;
t652 = t583 * t676;
t650 = t510 * t663;
t649 = t512 * t662;
t648 = t514 * t661;
t647 = t510 * t666;
t646 = t512 * t665;
t645 = t514 * t664;
t644 = t678 * t685;
t643 = t605 * t656;
t642 = t608 * t656;
t641 = t677 * t684;
t640 = t606 * t654;
t639 = t609 * t654;
t638 = t676 * t683;
t637 = t607 * t652;
t636 = t610 * t652;
t629 = t671 * t715 - t709;
t628 = t672 * t715 + t708;
t627 = t669 * t714 - t707;
t626 = t670 * t714 + t706;
t625 = t667 * t713 - t705;
t624 = t668 * t713 + t704;
t571 = pkin(2) * t586 + pkin(3) * t580;
t570 = pkin(2) * t585 + pkin(3) * t579;
t569 = pkin(2) * t584 + pkin(3) * t578;
t559 = t716 * t713;
t556 = t717 * t714;
t553 = t718 * t715;
t1 = [0, (t515 * t666 + t516 * t665 + t517 * t664) * t616, 0, 0, (t509 * t666 + t511 * t665 + t513 * t664 + (-t710 * t719 - t711 * t720 - t712 * t721) * t615) * t616, t608 * t647 + t609 * t646 + t610 * t645 + ((-t713 * t637 + (-t559 * t668 - t704 * t716) * t574) * t598 + (-t714 * t640 + (-t556 * t670 - t706 * t717) * t573) * t597 + (-t715 * t643 + (-t553 * t672 - t708 * t718) * t572) * t596) * t615, -t605 * t647 - t606 * t646 - t607 * t645 + ((-t713 * t636 + (-t559 * t667 + t705 * t716) * t574) * t598 + (-t714 * t639 + (-t556 * t669 + t707 * t717) * t573) * t597 + (-t715 * t642 + (-t553 * t671 + t709 * t718) * t572) * t596) * t615, 0; 0, (-t515 * t663 - t516 * t662 - t517 * t661) * t616, 0, 0, (-t509 * t663 - t511 * t662 - t513 * t661 + (t686 * t710 + t687 * t711 + t688 * t712) * t615) * t616, -t608 * t650 - t609 * t649 - t610 * t648 + ((t624 * t574 + t637 * t716) * t686 + (t626 * t573 + t640 * t717) * t687 + (t628 * t572 + t643 * t718) * t688) * t615, t605 * t650 + t606 * t649 + t607 * t648 + ((t625 * t574 + t636 * t716) * t686 + (t627 * t573 + t639 * t717) * t687 + (t629 * t572 + t642 * t718) * t688) * t615, 0; 0, (-t515 * t685 - t516 * t684 - t517 * t683) * t616, 0, 0, (-t509 * t685 - t511 * t684 - t513 * t683 + (t509 * t569 * t718 + t511 * t570 * t717 + t513 * t571 * t716) * t615) * t616, -t608 * t660 - t609 * t659 - t610 * t658 + ((t624 * t571 + t607 * t638) * t716 + (t626 * t570 + t606 * t641) * t717 + (t628 * t569 + t605 * t644) * t718) * t615, t605 * t660 + t606 * t659 + t607 * t658 + ((t625 * t571 + t610 * t638) * t716 + (t627 * t570 + t609 * t641) * t717 + (t629 * t569 + t608 * t644) * t718) * t615, 0;];
tau_reg  = t1;
