% Return the minimum parameter vector for
% P4PRRRR8V2G1A0
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% 
% Output:
% MPV [15x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = P4PRRRR8V2G1A0_minimal_parameter_para(pkin, m, mrSges, Ifges, koppelP)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6),zeros(4,3)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_minimal_parameter_para: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G1A0_minimal_parameter_para: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G1A0_minimal_parameter_para: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V2G1A0_minimal_parameter_para: Ifges has to be [4x6] (double)'); 
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_minimal_parameter_para: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From minimal_parameter_parrob_matlab.m
t1 = [m(1) + m(2) + m(3); Ifges(2,3) + Ifges(3,2) + 2 * pkin(6) * mrSges(3,3) + (pkin(2) ^ 2 + pkin(6) ^ 2) * m(3); m(3) * pkin(2) + mrSges(2,1); -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2); Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3); mrSges(3,1); mrSges(3,2); Ifges(4,3); mrSges(4,1); mrSges(4,2); m(4);];
MPV  = t1;
