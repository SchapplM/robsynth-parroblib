% Calculate Gravitation load for parallel robot
% P3PRRR1G2A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:19
% EndTime: 2020-03-09 21:18:19
% DurationCPUTime: 0.26s
% Computational Cost: add. (483->73), mult. (480->118), div. (45->5), fcn. (309->18), ass. (0->63)
t605 = pkin(7) + qJ(2,3);
t596 = qJ(3,3) + t605;
t584 = sin(t596);
t610 = legFrame(3,2);
t599 = sin(t610);
t602 = cos(t610);
t581 = t602 * g(1) - t599 * g(2);
t587 = cos(t596);
t613 = rSges(3,2) * g(3);
t615 = rSges(3,1) * g(3);
t557 = -t587 * (rSges(3,1) * t581 - t613) + (rSges(3,2) * t581 + t615) * t584;
t590 = sin(t605);
t593 = cos(t605);
t609 = g(3) * pkin(2) * m(3);
t614 = rSges(2,2) * g(3);
t616 = rSges(2,1) * g(3);
t630 = pkin(2) * t593;
t634 = 0.1e1 / (-t593 * t584 + t587 * t590);
t624 = t634 * (m(2) * (-rSges(2,1) * t581 + t614) * t593 + t590 * (t609 + m(2) * (rSges(2,2) * t581 + t616)) + (-t581 * t630 + t557) * m(3));
t637 = t584 * t624;
t606 = pkin(7) + qJ(2,2);
t597 = qJ(3,2) + t606;
t585 = sin(t597);
t611 = legFrame(2,2);
t600 = sin(t611);
t603 = cos(t611);
t582 = t603 * g(1) - t600 * g(2);
t588 = cos(t597);
t558 = -t588 * (rSges(3,1) * t582 - t613) + (rSges(3,2) * t582 + t615) * t585;
t591 = sin(t606);
t594 = cos(t606);
t629 = pkin(2) * t594;
t633 = 0.1e1 / (-t594 * t585 + t588 * t591);
t623 = t633 * (m(2) * (-rSges(2,1) * t582 + t614) * t594 + t591 * (t609 + m(2) * (rSges(2,2) * t582 + t616)) + (-t582 * t629 + t558) * m(3));
t636 = t585 * t623;
t607 = pkin(7) + qJ(2,1);
t598 = qJ(3,1) + t607;
t586 = sin(t598);
t612 = legFrame(1,2);
t601 = sin(t612);
t604 = cos(t612);
t583 = t604 * g(1) - t601 * g(2);
t589 = cos(t598);
t559 = -t589 * (rSges(3,1) * t583 - t613) + (rSges(3,2) * t583 + t615) * t586;
t592 = sin(t607);
t595 = cos(t607);
t628 = pkin(2) * t595;
t632 = 0.1e1 / (-t595 * t586 + t589 * t592);
t622 = t632 * (m(2) * (-rSges(2,1) * t583 + t614) * t595 + t592 * (t609 + m(2) * (rSges(2,2) * t583 + t616)) + (-t583 * t628 + t559) * m(3));
t635 = t586 * t622;
t631 = m(3) / pkin(3);
t627 = t557 * t634;
t626 = t558 * t633;
t625 = t559 * t632;
t621 = (pkin(2) * t590 + pkin(3) * t584) * t627;
t620 = (pkin(2) * t591 + pkin(3) * t585) * t626;
t619 = (pkin(2) * t592 + pkin(3) * t586) * t625;
t618 = 0.1e1 / pkin(2);
t608 = m(1) + m(2) + m(3);
t580 = t601 * g(1) + t604 * g(2);
t579 = t600 * g(1) + t603 * g(2);
t578 = t599 * g(1) + t602 * g(2);
t1 = [-m(4) * g(1) + (-t599 * t578 - t600 * t579 - t601 * t580) * t608 + (-t602 * t637 - t603 * t636 - t604 * t635 + (t602 * t621 + t603 * t620 + t604 * t619) * t631) * t618; -m(4) * g(2) + (-t578 * t602 - t579 * t603 - t580 * t604) * t608 + (t599 * t637 + t600 * t636 + t601 * t635 + (-t599 * t621 - t600 * t620 - t601 * t619) * t631) * t618; -m(4) * g(3) + (-t587 * t624 - t588 * t623 - t589 * t622 + ((pkin(3) * t589 + t628) * t625 + (pkin(3) * t588 + t629) * t626 + (pkin(3) * t587 + t630) * t627) * t631) * t618;];
taugX  = t1;
