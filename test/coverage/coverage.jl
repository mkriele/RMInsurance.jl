using Coverage

cd(joinpath(@__DIR__, "..", "..")) do
    Coveralls.submit(Coveralls.process_folder())
end
